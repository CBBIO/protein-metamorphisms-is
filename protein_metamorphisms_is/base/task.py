import ast
import multiprocessing
import pickle
import threading
import time
import traceback
from abc import ABC, abstractmethod

import pika
from pika import PlainCredentials

from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager


class BaseTaskInitializer(ABC):
    def __init__(self, conf, session_required=True):
        self.logger = setup_logger(self.__class__.__name__)
        self.logger.info(f"Initializing {self.__class__.__name__}")
        self.conf = conf
        self.processes = []
        self.threads = []

        if session_required:
            self.session_init()

    @property
    def extractor_queue(self):
        return f'{self.__class__.__name__.lower()}_extractor'

    @property
    def inserter_queue(self):
        return f'{self.__class__.__name__.lower()}_inserter'

    def setup_rabbitmq(self):
        credentials = PlainCredentials('guest', 'guest')  # Customize with your credentials
        self.connection = pika.BlockingConnection(
            pika.ConnectionParameters(host='localhost', credentials=credentials))
        self.channel = self.connection.channel()
        self.channel.queue_declare(queue=self.extractor_queue)
        self.channel.queue_declare(queue=self.inserter_queue)

    def start(self):
        try:
            self.setup_rabbitmq()
            self.enqueue()
            self.start_workers()
        finally:
            self.close_rabbitmq()

    def session_init(self):
        """
        Initialize the database session using DatabaseManager.

        Sets up the database connection and session using the DatabaseManager class.
        """
        self.logger.info("Initializing database session using DatabaseManager")
        db_manager = DatabaseManager(self.conf)
        self.engine = db_manager.get_engine()
        self.session = db_manager.get_session()

    def start(self):
        try:
            self.setup_rabbitmq()
            print('enq')
            self.enqueue()
            self.start_workers()
        finally:
            self.close_rabbitmq()

    def close_rabbitmq(self):
        try:
            if hasattr(self, 'channel') and self.channel.is_open:
                self.channel.close()
            if hasattr(self, 'connection') and self.connection.is_open:
                self.connection.close()
        except Exception as e:
            self.logger.error(f"Error closing RabbitMQ connection: {e}")

    def start_workers(self):
        # Start the processing workers
        for _ in range(self.conf['max_workers']):
            p = multiprocessing.Process(target=self.run_processor_worker)
            p.start()
            self.processes.append(p)

        # Start the DB insertion worker
        p = multiprocessing.Process(target=self.run_db_inserter_worker)
        p.start()
        self.processes.append(p)

        # Start the monitoring thread
        monitor_thread = threading.Thread(target=self.monitor_queues)
        monitor_thread.start()
        self.threads.append(monitor_thread)

        for p in self.processes:
            p.join()

        # Signal the monitor thread to stop
        self.stop_event.set()
        for t in self.threads:
            t.join()

    def run_processor_worker(self):
        self.setup_rabbitmq()
        self.channel.basic_qos(prefetch_count=1)
        self.channel.basic_consume(queue=self.extractor_queue, on_message_callback=self.callback)
        self.logger.info(f"Starting message consumption in processor worker for queue {self.extractor_queue}.")
        while not self.stop_event.is_set():
            try:
                self.connection.process_data_events(time_limit=1)
            except pika.exceptions.ConnectionClosedByBroker:
                break
            except pika.exceptions.AMQPChannelError as err:
                self.logger.error(f"Channel error: {err}")
                break
            except pika.exceptions.AMQPConnectionError:
                self.logger.error("Connection was closed, retrying...")
                break
        self.close_rabbitmq()

    def run_db_inserter_worker(self):
        self.setup_rabbitmq()
        self.channel.basic_qos(prefetch_count=1)
        self.channel.basic_consume(queue=self.inserter_queue, on_message_callback=self.db_inserter_callback)
        print(f"Starting message consumption in DB inserter worker for queue {self.inserter_queue}.")
        while not self.stop_event.is_set():
            try:
                self.connection.process_data_events(time_limit=1)
            except pika.exceptions.ConnectionClosedByBroker:
                break
            except pika.exceptions.AMQPChannelError as err:
                self.logger.error(f"Channel error: {err}")
                break
            except pika.exceptions.AMQPConnectionError:
                self.logger.error("Connection was closed, retrying...")
                break
        self.close_rabbitmq()

    def monitor_queues(self):
        self.setup_rabbitmq()
        while not self.stop_event.is_set():
            queues = [self.extractor_queue, self.inserter_queue]
            empty_queues = 0
            for queue in queues:
                queue_state = self.channel.queue_declare(queue=queue, passive=True)
                if queue_state.method.message_count == 0:
                    empty_queues += 1
            if empty_queues == len(queues):
                self.stop_event.set()
                break
            time.sleep(5)
        self.close_rabbitmq()

    def publish_task(self, data):
        if not isinstance(data, bytes):
            data = str(data).encode('utf-8')
        self.channel.basic_publish(exchange='', routing_key=self.extractor_queue, body=data)

    def db_inserter_callback(self, ch, method, properties, body):
        try:
            # Deserializar el registro con pickle
            record = pickle.loads(body)
            self.store_entry(record)
            ch.basic_ack(delivery_tag=method.delivery_tag)
        except Exception as e:
            error_message = f"Failed to insert record into DB: {e}"
            self.logger.error(error_message)
            ch.basic_nack(delivery_tag=method.delivery_tag, requeue=False)

    def callback(self, ch, method, properties, body):
        body_decoded = body.decode('utf-8')

        # Convierte el string de JSON a un diccionario
        data = ast.literal_eval(body_decoded)
        reference = data[self.reference_attribute]
        try:
            self.logger.info(f"Processing message for {self.reference_attribute}: {reference}")
            record = self.process(reference)
            if record:
                # Serialize the record with pickle
                pickled_record = pickle.dumps(record)
                # Publish the serialized record to the DB insertion queue
                self.channel.basic_publish(exchange='', routing_key=self.inserter_queue, body=pickled_record)
                self.logger.debug(
                    f"Published pickled record to inserter queue for {self.reference_attribute} code: {reference}")
            ch.basic_ack(delivery_tag=method.delivery_tag)
        except ValueError as e:
            self.logger.warning(f"No record found for {self.reference_attribute} code {reference}: {e}")
            ch.basic_ack(delivery_tag=method.delivery_tag)  # Acknowledge the message to remove it from the queue
        except Exception as e:
            self.logger.error(
                f"Failed to process message for {self.reference_attribute} code {reference}: {e}\n{traceback.format_exc()}")
            ch.basic_nack(delivery_tag=method.delivery_tag, requeue=False)

    @abstractmethod
    def process(self, target):
        pass

    @abstractmethod
    def enqueue(self):
        pass

    @abstractmethod
    def store_entry(self, record):
        pass
