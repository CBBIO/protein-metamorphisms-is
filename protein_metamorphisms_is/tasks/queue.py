"""
Queue Tasks
===========

Queue tasks handle complex and large-scale processes by leveraging RabbitMQ as a messaging system.
These tasks are designed to be distributed across multiple workers, ensuring efficient task processing
and fault tolerance.

**Purpose**

The `QueueTaskInitializer` class extends the `BaseTaskInitializer` and provides the necessary framework for
 managing tasks distributed via RabbitMQ. It initializes the RabbitMQ connection, manages worker processes,
 and ensures that tasks are processed in a reliable, scalable, and distributed manner.

**Customization**

To create a custom queue-based task, subclass `QueueTaskInitializer` and implement the `enqueue`, `process`, and
 `store_entry` methods. These methods define the logic for enqueuing tasks, processing them, and storing the results
 in the database.

**Key Features**

- **RabbitMQ Integration**: Handles the setup and management of RabbitMQ queues for task distribution.
- **Worker Management**: Manages multiple worker processes for parallel task execution.
- **Logger Inheritance**: Inherits task-specific logging from `BaseTaskInitializer`.
- **Database Integration**: Inherits database session management from `BaseTaskInitializer`.
- **Extensibility**: Abstract methods are provided to allow developers to define custom task logic.

**Example Usage**

Here is an example of how to subclass `QueueTaskInitializer`:

.. code-block:: python

   from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer

   class MyCustomQueueTask(QueueTaskInitializer):
       def enqueue(self):
           # Implementation of the task enqueuing logic
           pass

       def process(self, target):
           # Processing logic for the target data
           pass

       def store_entry(self, record):
           # Logic to store the processed record in the database
           pass
"""

import multiprocessing
import pickle
import threading
import time
from abc import abstractmethod
import pika
import requests
from pika import PlainCredentials

from protein_metamorphisms_is.tasks.base import BaseTaskInitializer


class QueueTaskInitializer(BaseTaskInitializer):
    """
    The QueueTaskInitializer class extends BaseTaskInitializer to manage tasks
    that are distributed and processed via RabbitMQ queues.

    It offers the necessary infrastructure to configure RabbitMQ, coordinate
    worker processes, and ensure that task6s are handled efficiently and reliably.

    Attributes:
        processes (list): A list of multiprocessing.Process objects representing the worker processes.
        threads (list): A list of threading.Thread objects for auxiliary tasks like monitoring queues.
        stop_event (multiprocessing.Event): An event to signal workers and threads to stop.
        connection_params (pika.ConnectionParameters): Parameters for connecting to the RabbitMQ server.
    """

    def __init__(self, conf, session_required=True):
        """
        Initialize the QueueTaskInitializer.

        This constructor initializes the configuration, and if required, sets up
        a database session. It also prepares the RabbitMQ connection parameters and
        initializes the processes and threads lists.

        Args:
            conf (dict): Configuration dictionary.
            session_required (bool): Whether a database session is required.
                                     If True, the session is initialized.
        """
        super().__init__(conf, session_required)
        self.processes = []
        self.threads = []
        self.stop_event = multiprocessing.Event()
        self.connection_params = pika.ConnectionParameters(
            host=self.conf['rabbitmq_host'],
            heartbeat=300,
            credentials=PlainCredentials(self.conf['rabbitmq_user'], self.conf['rabbitmq_password'])
        )

    def setup_rabbitmq(self):
        """
        Set up RabbitMQ by declaring the necessary queues.

        This method connects to RabbitMQ using the provided credentials and declares
        the necessary queues for task computing and data insertion.

        Raises:
            Exception: If there is an issue setting up RabbitMQ.
        """
        try:
            credentials = PlainCredentials(self.conf['rabbitmq_user'], self.conf['rabbitmq_password'])
            connection = pika.BlockingConnection(
                pika.ConnectionParameters(host=self.conf['rabbitmq_host'], credentials=credentials))
            channel = connection.channel()

            # Declare necessary queues
            channel.queue_declare(queue=self.computing_queue, durable=True)
            channel.queue_declare(queue=self.inserting_queue, durable=True)

        except Exception as e:
            raise e

    @property
    def computing_queue(self):
        """
        Property to get the name of the extractor queue.

        The extractor queue is used for processing tasks that extract data.

        Returns:
            str: Name of the extractor queue.
        """
        return f'{self.__class__.__name__.lower()}_computing'

    @property
    def inserting_queue(self):
        """
        Property to get the name of the inserter queue.

        The inserter queue is used for tasks that insert processed data into the database.

        Returns:
            str: Name of the inserter queue.
        """
        return f'{self.__class__.__name__.lower()}_inserting'

    def start(self):
        """
        Start the task processing pipeline.

        This method sets up RabbitMQ, enqueues tasks, and starts worker processes.
        It also handles the graceful shutdown of workers when the stop event is set.
        """
        try:
            self.setup_rabbitmq()
            self.enqueue()
            self.start_workers()
        finally:
            self.stop_event.set()

    def start_workers(self):
        """
        Start the worker processes for task computing and database insertion of the results.

        This method spawns multiple worker processes for processing tasks and
        inserts processed data into the database. It also starts a monitoring thread
        to oversee the queues.
        """
        try:
            for _ in range(self.conf['max_workers']):
                p = multiprocessing.Process(target=self.run_processor_worker, args=(self.stop_event,))
                p.start()
                self.processes.append(p)

            p = multiprocessing.Process(target=self.run_db_inserter_worker, args=(self.stop_event,))
            p.start()
            self.processes.append(p)

            monitor_thread = threading.Thread(target=self.monitor_queues)
            monitor_thread.start()
            self.threads.append(monitor_thread)

            for p in self.processes:
                p.join()
        finally:
            self.stop_event.set()
            for t in self.threads:
                t.join()

    def run_processor_worker(self, stop_event):
        """
        Run the processor worker.

        This worker consumes messages from the computing_queue and processes
        the tasks. It acknowledges messages after processing.

        Args:
            stop_event (multiprocessing.Event): Event to signal the worker to stop.
        """
        with self.create_rabbitmq_connection() as connection:
            channel = connection.channel()
            channel.basic_qos(prefetch_count=1)
            channel.basic_consume(queue=self.computing_queue, on_message_callback=self.callback)
            self.consume_messages(channel, stop_event)

    def run_db_inserter_worker(self, stop_event):
        """
        Run the database inserter worker.

        This worker consumes messages from the inserting_queue and save
        the processed data into the database.

        Args:
            stop_event (multiprocessing.Event): Event to signal the worker to stop.
        """
        with self.create_rabbitmq_connection() as connection:
            channel = connection.channel()
            channel.basic_qos(prefetch_count=1)
            channel.basic_consume(queue=self.inserting_queue, on_message_callback=self.db_inserter_callback)
            self.consume_messages(channel, stop_event)

    def create_rabbitmq_connection(self):
        """
        Create a connection to RabbitMQ.

        This method establishes a connection to RabbitMQ using the provided
        connection parameters.

        Returns:
            pika.BlockingConnection: The RabbitMQ connection.
        """
        return pika.BlockingConnection(self.connection_params)

    def consume_messages(self, channel, stop_event):
        """
        Consume messages from a RabbitMQ queue.

        This method processes messages from the given channel and acknowledges them.
        It stops consuming when the stop event is set.

        Args:
            channel (pika.channel.Channel): The RabbitMQ channel to consume from.
            stop_event (multiprocessing.Event): Event to signal when to stop consuming messages.
        """
        while not stop_event.is_set():
            try:
                channel.connection.process_data_events(time_limit=1)
            except pika.exceptions.ConnectionClosedByBroker:
                break
            except pika.exceptions.AMQPChannelError:
                break
            except pika.exceptions.AMQPConnectionError:
                break

    def check_messages_in_memory(self, queue_name):
        """
        Check the number of messages in memory for a given queue.

        This method accesses the RabbitMQ management API to retrieve the number of messages
        currently stored in memory for the specified queue.

        Args:
            queue_name (str): The name of the queue to check.

        Returns:
            int: The number of messages in memory for the queue.
        """
        host = self.conf['rabbitmq_host']
        user = self.conf['rabbitmq_user']
        password = self.conf['rabbitmq_password']
        url = f'http://{host}:15672/api/queues/%2F/{queue_name}'
        auth = requests.auth.HTTPBasicAuth(user, password)

        try:
            response = requests.get(url, auth=auth)
            response.raise_for_status()

            queue_info = response.json()
            messages_in_memory = queue_info.get('messages_ram', 0)

            if isinstance(messages_in_memory, str):
                messages_in_memory = int(messages_in_memory)

            return messages_in_memory

        except requests.exceptions.RequestException:
            return 0

        except ValueError:
            return 0

    def monitor_queues(self):
        """
        Monitor the RabbitMQ queues.

        This method periodically checks the status of the queues to determine
        if they are empty. If all queues are empty, it stops the workers.

        This method runs in a separate thread and monitors the queues every 5 seconds.

        """
        while not self.stop_event.is_set():
            time.sleep(5)

            try:
                with self.create_rabbitmq_connection() as connection:
                    channel = connection.channel()
                    queues = [self.computing_queue, self.inserting_queue]
                    all_queues_empty = True
                    all_memory_queues_empty = True

                    for queue in queues:
                        try:
                            queue_state = channel.queue_declare(queue=queue, passive=True)

                            if queue_state.method.message_count > 0:
                                all_queues_empty = False

                            messages_in_memory = self.check_messages_in_memory(queue)
                            if messages_in_memory > 0:
                                all_memory_queues_empty = False

                        except pika.exceptions.ChannelClosedByBroker:
                            self.stop_event.set()
                            return

                    if all_queues_empty and all_memory_queues_empty:
                        self.stop_event.set()
                        break
                    else:
                        continue

            except pika.exceptions.AMQPConnectionError:
                self.stop_event.set()
                break
            except Exception:
                self.stop_event.set()
                break

    def publish_task(self, data):
        """
        Publish a task to the extractor queue.

        This method serializes the data processing and publishes it to the computing queue.

        Args:
            data (any): The task data to be published.
        """
        if not isinstance(data, bytes):
            data = pickle.dumps(data)
        with self.create_rabbitmq_connection() as connection:
            channel = connection.channel()
            channel.basic_publish(exchange='', routing_key=self.computing_queue, body=data)

    def db_inserter_callback(self, ch, method, properties, body):
        """
        Callback function for the database inserter worker.

        This function is triggered when a message is received in the inserter queue.
        It deserializes the message, stores the entry in the database, and acknowledges the message.

        Args:
            ch (pika.channel.Channel): The channel object.
            method (pika.spec.Basic.Deliver): The method frame with delivery details.
            properties (pika.spec.BasicProperties): The properties of the message.
            body (bytes): The body of the message.
        """
        try:
            record = pickle.loads(body)
            self.store_entry(record)
            ch.basic_ack(delivery_tag=method.delivery_tag)
        except Exception:
            ch.basic_nack(delivery_tag=method.delivery_tag, requeue=False)

    def callback(self, ch, method, properties, body):
        """
        Callback function for the processor worker.

        This function is triggered when a message is received in the extractor queue.
        It processes the message, and if a record is generated, it publishes the record
        to the inserter queue.

        Args:
            ch (pika.channel.Channel): The channel object.
            method (pika.spec.Basic.Deliver): The method frame with delivery details.
            properties (pika.spec.BasicProperties): The properties of the message.
            body (bytes): The body of the message.
        """
        try:
            data = pickle.loads(body)
            reference = data
            record = self.process(reference)
            if record:
                record_bytes = pickle.dumps(record)
                ch.basic_publish(exchange='', routing_key=self.inserting_queue, body=record_bytes)

            ch.basic_ack(delivery_tag=method.delivery_tag)
        except ValueError:
            ch.basic_ack(delivery_tag=method.delivery_tag)
        except Exception:
            ch.basic_nack(delivery_tag=method.delivery_tag, requeue=True)

    @abstractmethod
    def enqueue(self):
        """
        Abstract method to enqueue tasks.

        This method should be implemented to retrieve from the database the data source and serialize it to the
        computing queue.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass

    @abstractmethod
    def process(self, target):
        """
        Abstract method to process tasks.

        This method should be implemented by subclasses to define how each task
        is processed.

        Args:
            target: The target data to be processed.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass

    @abstractmethod
    def store_entry(self, record):
        """
        Abstract method to store processed entries.

        This method should be implemented by subclasses to define how processed
        records are stored in the database.

        Args:
            record: The processed data record to be stored.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass
