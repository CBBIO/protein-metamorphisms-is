"""
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

- **Session Management**: Integrates seamlessly with the database to manage sessions
  and maintain data consistency.
- **Configuration Handling**: Loads and processes configuration constants from YAML files
  to ensure that all tasks are initialized with the correct settings.
- **Extensibility**: Abstract methods are provided to allow developers to define
  specific task logic for their bioinformatics workflows.


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
"""

from abc import ABC, abstractmethod
from pathlib import Path

import yaml
from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager
from protein_metamorphisms_is.sql.constants import handle_sequence_embedding_types, handle_structural_alignment_types


class BaseTaskInitializer(ABC):
    """
    The BaseTaskInitializer class provides a foundation for creating
    bioinformatics tasks that interact with the database via ORM.

    This class is abstract and should be subclassed to define specific task
    processing logic. It handles session management, loading configuration constants,
    and provides an interface for starting, processing, and storing task data.

    Additionally, it initializes a task-specific logger.

    Attributes:
        conf (dict): Configuration dictionary loaded from YAML or other sources.
        logger (Logger): Logger instance for logging task-specific information.
        session (Session): Database session used for ORM operations.
    """

    def __init__(self, conf, session_required=True):
        """
        Initialize the BaseTaskInitializer.

        This constructor sets up the logger for the task, initializes the configuration,
        and, if required, sets up a database session and loads constants.

        Args:
            conf (dict): Configuration dictionary.
            session_required (bool): Whether a database session is required.
                                     If True, the session is initialized.

        Raises:
            Exception: If there is an issue initializing the session or loading constants.
        """
        self.logger = setup_logger(self.__class__.__name__, conf['log_path'])
        self.logger.info(f"Initializing {self.__class__.__name__}")
        self.conf = conf

        if session_required:
            self.session_init()
            self.load_constants(conf['constants'])

    def session_init(self):
        """
        Initialize the database session using DatabaseManager.

        Sets up the database connection and session using the DatabaseManager class.
        Logs the initialization process.

        Raises:
            Exception: If the session initialization fails.
        """
        self.logger.info("Initializing database session using DatabaseManager")
        db_manager = DatabaseManager(self.conf)
        self.engine = db_manager.get_engine()
        self.session = db_manager.get_session()
        self.pool = db_manager.get_pool()

    def load_constants(self, constants_path):
        """
        Load and handle predefined constants.

        This method processes the constants YAML file and ensures that all
        required structural and prediction constants are up-to-date in the database.

        Logs the process of loading constants.

        Args:
            constants_path (str): Path to the YAML file containing configuration constants.

        Raises:
            FileNotFoundError: If the constants file cannot be found.
            yaml.YAMLError: If there is an error parsing the YAML file.
        """
        try:
            self.logger.info(f"Loading constants from {constants_path}")

            # Verify file exists
            if not Path(constants_path).is_file():
                raise FileNotFoundError(f"Constants file not found at: {constants_path}")

            # Load and parse YAML safely
            with open(constants_path, 'r', encoding='utf-8') as f:
                constants = yaml.safe_load(f)

            # Handle specific constant types
            handle_structural_alignment_types(self.session, constants)
            handle_sequence_embedding_types(self.session, constants)

        except yaml.YAMLError as e:
            self.logger.error(f"Error parsing YAML file: {str(e)}")
            raise
        except Exception as e:
            self.logger.error(f"Error loading constants: {str(e)}")
            raise

    @abstractmethod
    def start(self):
        """
        Start the operation over the data process.

        This method should be implemented by all subclasses to define
        the specific data operation logic for each bioinformatics data source.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass

    @abstractmethod
    def process(self, target):
        """
        Process the given target.

        This method should be implemented by all subclasses to define
        the specific processing logic for each bioinformatics data source.

        Args:
            target: The target data to be processed.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass

    @abstractmethod
    def store_entry(self, record):
        """
        Store the processed entry.

        This method should be implemented by all subclasses to define
        the specific storage logic for processed data.

        Args:
            record: The processed data record to be stored.

        Raises:
            NotImplementedError: If the subclass does not implement this method.
        """
        pass
