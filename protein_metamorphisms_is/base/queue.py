import ast
import multiprocessing
import json
import threading
import time
import traceback
from abc import abstractmethod

import pika
from pika import PlainCredentials
from retry import retry

from protein_metamorphisms_is.base.base import BaseTaskInitializer
from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager

class QueueTaskInitializer(BaseTaskInitializer):
    def __init__(self, conf, session_required=True):
        super().__init__(conf, session_required)
        self.processes = []
        self.threads = []
        self.stop_event = multiprocessing.Event()
        self.connection = None
        self.channel = None

    @property
    def extractor_queue(self):
        return f'{self.__class__.__name__.lower()}_extractor'

    @property
    def inserter_queue(self):
        return f'{self.__class__.__name__.lower()}_inserter'

    def setup_rabbitmq(self):
        try:
            credentials = PlainCredentials(self.conf['rabbitmq_user'], self.conf['rabbitmq_password'])
            self.connection = pika.BlockingConnection(
                pika.ConnectionParameters(host=self.conf['rabbitmq_host'], credentials=credentials))
            self.channel = self.connection.channel()
            self.channel.queue_declare(queue=self.extractor_queue)
            self.channel.queue_declare(queue=self.inserter_queue)
        except Exception as e:
            self.logger.error(f"Error setting up RabbitMQ: {e}")
            self.close_rabbitmq()
            raise

    def start(self):
        try:
            self.setup_rabbitmq()
            self.enqueue()
            self.start_workers()
        finally:
            self.close_rabbitmq()

    def close_rabbitmq(self):
        try:
            if self.channel and self.channel.is_open:
                self.channel.close()
            if self.connection and self.connection.is_open:
                self.connection.close()
        except Exception as e:
            self.logger.error(f"Error closing RabbitMQ connection: {e}")

    def start_workers(self):
        try:
            for _ in range(self.conf['max_workers']):
                p = multiprocessing.Process(target=self.run_processor_worker, args=(self.stop_event,))
                p.start()
                self.processes.append(p)

            p = multiprocessing.Process(target=self.run_db_inserter_worker, args=(self.stop_event,))
            p.start()
            self.processes.append(p)

            # monitor_thread = threading.Thread(target=self.monitor_queues)
            # monitor_thread.start()
            # self.threads.append(monitor_thread)

            for p in self.processes:
                p.join()
        except Exception as e:
            self.logger.error(f"Error starting workers: {e}")
        finally:
            self.stop_event.set()
            for t in self.threads:
                t.join()

    def run_processor_worker(self, stop_event):
        with self._create_rabbitmq_connection() as channel:
            channel.basic_qos(prefetch_count=1)
            channel.basic_consume(queue=self.extractor_queue, on_message_callback=self.callback)
            self.logger.info(f"Starting message consumption in processor worker for queue {self.extractor_queue}.")
            self._consume_messages(channel, stop_event)

    def run_db_inserter_worker(self, stop_event):
        with self._create_rabbitmq_connection() as channel:
            channel.basic_qos(prefetch_count=1)
            channel.basic_consume(queue=self.inserter_queue, on_message_callback=self.db_inserter_callback)
            self.logger.info(f"Starting message consumption in DB inserter worker for queue {self.inserter_queue}.")
            self._consume_messages(channel, stop_event)

    def _create_rabbitmq_connection(self):
        credentials = PlainCredentials(self.conf['rabbitmq_user'], self.conf['rabbitmq_password'])
        connection = pika.BlockingConnection(
            pika.ConnectionParameters(host=self.conf['rabbitmq_host'], credentials=credentials))
        channel = connection.channel()
        return channel

    def _consume_messages(self, channel, stop_event):
        while not stop_event.is_set():
            try:
                channel.connection.process_data_events(time_limit=1)
            except pika.exceptions.ConnectionClosedByBroker:
                self.logger.error("Connection closed by broker, stopping worker.")
                break
            except pika.exceptions.AMQPChannelError as err:
                self.logger.error(f"Channel error: {err}")
                break
            except pika.exceptions.AMQPConnectionError:
                self.logger.error("Connection error, retrying...")
                break

    def monitor_queues(self):
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

        self.logger.info("Queues are empty, stopping all workers.")
        self.stop_event.set()

    def publish_task(self, data):
        if not isinstance(data, str):
            data = json.dumps(data)
        if not self.channel or not self.channel.is_open:
            self.setup_rabbitmq()
        self.channel.basic_publish(exchange='', routing_key=self.extractor_queue, body=data)

    def db_inserter_callback(self, ch, method, properties, body):
        self.logger.info("Message received in DB inserter worker")
        try:
            record = json.loads(body)
            self.store_entry(record)
            ch.basic_ack(delivery_tag=method.delivery_tag)
        except Exception as e:
            error_message = f"Failed to insert record into DB: {e}"
            self.logger.error(error_message)
            ch.basic_nack(delivery_tag=method.delivery_tag, requeue=False)

    def callback(self, ch, method, properties, body):
        body_decoded = body.decode('utf-8')
        data = json.loads(body_decoded)
        reference = data
        try:
            self.logger.info(f"Processing message for {self.reference_attribute}: {reference}")
            record = self.process(reference)
            if record:
                record_json = json.dumps(record)
                if not self.channel or not self.channel.is_open:
                    self.setup_rabbitmq()
                self.channel.basic_publish(exchange='', routing_key=self.inserter_queue, body=record_json)
                self.logger.debug(
                    f"Published record to inserter queue for {self.reference_attribute} code: {reference}")
            ch.basic_ack(delivery_tag=method.delivery_tag)
        except ValueError as e:
            self.logger.warning(f"No record found for {self.reference_attribute} code {reference}: {e}")
            ch.basic_ack(delivery_tag=method.delivery_tag)
        except Exception as e:
            self.logger.error(
                f"Failed to process message for {self.reference_attribute} code {reference}: {e}\n{traceback.format_exc()}")
            ch.basic_nack(delivery_tag=method.delivery_tag, requeue=False)

    @abstractmethod
    def enqueue(self):
        pass

    @abstractmethod
    def process(self, target):
        pass

    @abstractmethod
    def store_entry(self, record):
        pass
