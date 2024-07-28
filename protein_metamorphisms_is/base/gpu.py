import importlib
import multiprocessing
import pickle
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
        last_model_type = None  # Track the last model type loaded
        with self._create_rabbitmq_connection() as channel:
            queue_name = f"{self.extractor_queue}_{model_type}"
            channel.basic_qos(prefetch_count=1)
            while not self.stop_event.is_set():
                method_frame, header_frame, body = channel.basic_get(queue=queue_name, auto_ack=False)
                if method_frame:
                    # Check if model needs to be changed
                    if model_type != last_model_type:
                        if last_model_type is not None:
                            self.unload_model(last_model_type)  # Unload previous model
                        self.load_model(model_type)  # Load the required model
                        last_model_type = model_type

                    # Process task
                    self.callback(channel, method_frame, header_frame, body)
                else:
                    break
            if last_model_type is not None:
                self.unload_model(last_model_type)  # Ensure last model is unloaded

        self.logger.info(f"Finished processing for queue {queue_name}.")

    def load_model(self, model_type):
        """Load the model into memory."""
        type_obj = self.types[model_type]
        module = type_obj['module']
        model = module.load_model(type_obj['model_name'])
        tokenizer = module.load_tokenizer(type_obj['model_name'])
        self.model_instances[model_type] = model
        self.tokenizer_instances[model_type] = tokenizer

    def unload_model(self, model_type):
        """Unload the model from memory."""
        if model_type in self.model_instances:
            del self.model_instances[model_type]
            del self.tokenizer_instances[model_type]



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
        if not isinstance(batch_data, bytes):  # Asegúrate de que los datos sean bytes antes de publicarlos
            batch_data = pickle.dumps(batch_data)
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
