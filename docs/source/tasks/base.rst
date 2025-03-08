Base Tasks
==========

Base tasks form the core of the system’s operations, providing a foundational framework for data extraction,
processing, and storage. These tasks can be implemented in both single-threaded and multi-threaded or
multiprocessing environments. They are directly integrated with the database via Object-Relational Mapping (ORM),
facilitating session management and ensuring data consistency.

**Purpose**

The `BaseTaskInitializer` class serves as an abstract base class that defines the common structure and behavior
for all base tasks within the system. It provides essential methods for initializing database sessions,
handling configuration constants, and defining abstract methods that must be implemented by subclasses to process specific bioinformatics data.

**Customization**

To create a custom task, subclass `BaseTaskInitializer` and implement the `start`, `process`, and `store_entry`
methods. These methods define the logic for processing specific data sources and storing the processed data
in the database.

**Key Features**

- **Session Management**: Integrates seamlessly with the database to manage sessions and maintain data consistency.
- **Configuration Handling**: Loads and processes configuration constants from YAML files to ensure that all tasks are initialized with the correct settings.
- **Extensibility**: Abstract methods are provided to allow developers to define specific task logic for their bioinformatics workflows.

**Example Usage**

Here is an example of how to subclass `BaseTaskInitializer`:

.. code-block:: python

   from protein_metamorphisms_is.tasks.base import BaseTaskInitializer

   class MyCustomTask(BaseTaskInitializer):
       def start(self):
           # Implementation of the start method
           pass

       def process(self, target):
           # Processing logic for the target data
           pass

       def store_entry(self, record):
           # Logic to store the processed record in the database
           pass