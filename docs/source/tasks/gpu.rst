GPU Tasks
=========

GPU tasks are specialized for managing and executing computationally intensive operations using Graphics
Processing Units (GPUs). These tasks are optimized to efficiently handle large-scale data processing,
particularly during model training and inference, where high-performance computing is critical.

**Purpose**

The `GPUTaskInitializer` class extends the `QueueTaskInitializer` to provide the necessary
framework for managing GPU-accelerated tasks. It handles the dynamic loading and unloading of GPU models,
manages worker processes, and ensures that tasks are processed efficiently using GPU resources.

**Customization**

To create a custom GPU-based task, subclass `GPUTaskInitializer` and implement the `enqueue`, `process`, and
 `store_entry` methods. These methods define the logic for enqueuing tasks, processing them using GPUs, and storing the results in the database.

**Key Features**

- **Efficient Resource Management**: Dynamically load and unload GPU models to maximize resource utilization.
- **Scalable Task Execution**: Handle large batches of data, optimizing throughput during model training and inference.
- **Seamless Integration**: Fully integrates with RabbitMQ for distributed task management, ensuring reliable and scalable task execution.
- **Worker Management**: Manages GPU-specific worker processes to handle tasks in parallel.
- **Extensibility**: Abstract methods are provided to allow developers to define custom GPU task logic.

**Example Usage**

Here is an example of how to subclass `GPUTaskInitializer`:

.. code-block:: python

   from protein_metamorphisms_is.tasks.gpu import GPUTaskInitializer

   class MyCustomGPUTask(GPUTaskInitializer):
       def enqueue(self):
           # Implementation of the task enqueuing logic for GPU processing
           pass

       def process(self, target):
           # Processing logic for the target data using GPU
           pass

       def store_entry(self, record):
           # Logic to store the processed record in the database
           pass