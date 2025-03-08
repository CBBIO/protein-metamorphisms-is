Queue Tasks
=============

Queue tasks extend the core functionality of the system by managing distributed task processing
through RabbitMQ queues. These tasks allow scalable, asynchronous processing of data, ensuring
reliable execution and optimized resource utilization.

**Purpose**

The QueueTaskInitializer class extends BaseTaskInitializer to handle task distribution via RabbitMQ.
It provides methods for setting up queues, managing worker processes, and handling message consumption
and processing.

**Customization**

To implement a queue-based task, subclass QueueTaskInitializer and define the enqueue, process,
and store_entry methods. These methods specify how tasks are enqueued, processed, and stored.

**Key Features**

- **RabbitMQ Integration**: Configures RabbitMQ connections and manages message queues efficiently.
- **Worker Management**: Supports both multiprocessing and multithreading for handling task execution.
- **Scalability**: Ensures efficient load balancing and queue monitoring for distributed systems.
- **Reliable Task Execution**: Uses acknowledgments and error handling to ensure robust message processing.
- **Extensibility**: Provides abstract methods to allow customization for various distributed processing needs.

**Example Usage**

Here is an example of how to subclass `QueueTaskInitializer`:

.. code-block:: python

    from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer

    class MyQueueTask(QueueTaskInitializer):
        def enqueue(self):
            # Task enqueueing logic
            pass

        def process(self, target):
            # Processing logic for the task
            pass

        def store_entry(self, record):
            # Logic to store processed records in the database
            pass