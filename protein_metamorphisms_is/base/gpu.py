import multiprocessing
import threading
import time
import json

import pika
from pika import PlainCredentials
from protein_metamorphisms_is.base.queue import QueueTaskInitializer

class GPUTaskInitializer(QueueTaskInitializer):
    def __init__(self, conf, session_required=True):
        super().__init__(conf, session_required)
        self.stop_event = multiprocessing.Event()

    def setup_rabbitmq(self):
        try:
            credentials = PlainCredentials(self.conf['rabbitmq_user'], self.conf['rabbitmq_password'])
            self.connection = pika.BlockingConnection(
                pika.ConnectionParameters(host=self.conf['rabbitmq_host'], credentials=credentials))
            self.channel = self.connection.channel()

            # Declare queues for each embedding type
            for model_type in self.conf['embedding']['types']:
                queue_name = f"{self.extractor_queue}_{model_type}"
                self.channel.queue_declare(queue=queue_name)
            self.channel.queue_declare(queue=self.inserter_queue)

        except Exception as e:
            self.logger.error(f"Error setting up RabbitMQ: {e}")
            self.close_rabbitmq()
            raise

    def start_workers(self):
        try:
            # Setup RabbitMQ and declare queues
            self.setup_rabbitmq()

            # Start the database inserter process
            db_inserter_process = multiprocessing.Process(target=self.run_db_inserter_worker)
            db_inserter_process.start()
            self.processes.append(db_inserter_process)

            # Start the monitor thread
            monitor_thread = threading.Thread(target=self.monitor_queues)
            monitor_thread.start()
            self.threads.append(monitor_thread)

            # Start the processing for each model type sequentially
            for model_type in self.conf['embedding']['types']:
                self.run_processor_worker_sequential(model_type)

            # Wait for the database inserter process to finish
            db_inserter_process.join()

        except Exception as e:
            self.logger.error(f"Error starting workers: {e}")

        finally:
            self.cleanup()

    def cleanup(self):
        self.stop_event.set()
        for t in self.threads:
            t.join()
        for p in self.processes:
            p.terminate()  # Ensure all processes are terminated properly

    def run_processor_worker_sequential(self, model_type):
        with self._create_rabbitmq_connection() as channel:
            queue_name = f"{self.extractor_queue}_{model_type}"
            channel.basic_qos(prefetch_count=1)
            while not self.stop_event.is_set():
                method_frame, header_frame, body = channel.basic_get(queue=queue_name, auto_ack=False)
                if method_frame:
                    self.callback(channel, method_frame, header_frame, body)
                else:
                    break
            self.logger.info(f"Finished processing for queue {queue_name}.")

    def run_db_inserter_worker(self):
        with self._create_rabbitmq_connection() as channel:
            channel.basic_qos(prefetch_count=1)
            channel.basic_consume(queue=self.inserter_queue, on_message_callback=self.db_inserter_callback)
            self.logger.info(f"Starting message consumption in DB inserter worker for queue {self.inserter_queue}.")
            self._consume_messages(channel, self.stop_event)

    def monitor_queues(self):
        with self._create_rabbitmq_connection() as channel:
            while not self.stop_event.is_set():
                queues = [f"{self.extractor_queue}_{model_type}" for model_type in self.conf['embedding']['types']]
                empty_queues = 0
                for queue in queues:
                    queue_state = channel.queue_declare(queue=queue, passive=True)
                    if queue_state.method.message_count == 0:
                        empty_queues += 1

                if empty_queues == len(queues):
                    self.stop_event.set()
                    break
                time.sleep(5)

            self.logger.info("Queues are empty, stopping all workers.")
            self.stop_event.set()

    def publish_task(self, batch_data, model_type):
        print(batch_data)
        if not isinstance(batch_data, str):
            batch_data = json.dumps(batch_data)
        if not self.channel or not self.channel.is_open:
            self.setup_rabbitmq()
        queue_name = f"{self.extractor_queue}_{model_type}"
        self.channel.basic_publish(exchange='', routing_key=queue_name, body=batch_data)

    def _create_rabbitmq_connection(self):
        credentials = PlainCredentials(self.conf['rabbitmq_user'], self.conf['rabbitmq_password'])
        connection = pika.BlockingConnection(
            pika.ConnectionParameters(host=self.conf['rabbitmq_host'], credentials=credentials))
        channel = connection.channel()
        return channel
