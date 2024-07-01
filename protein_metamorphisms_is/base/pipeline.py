import multiprocessing
import queue
import time
from abc import ABC, abstractmethod

from redis import Redis
from rq import Queue

from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager


class PipelineBase(ABC):
    """
    An abstract base class for operating bioinformatics data.

    This class provides a framework for connecting to and interacting with various
    bioinformatics data sources. It is designed to be subclassed with specific
    implementations for different data sources.

    Attributes:
        conf (dict): Configuration parameters for the extractor.
        logger (Logger): Logger object for logging information.
        session (Session, optional): A SQLAlchemy session for database interactions.

    Args:
        conf (dict): Configuration dictionary containing necessary parameters.
        session_required (bool): Flag to indicate if a database session is required.
    """

    def __init__(self, conf, session_required=False):
        """
        Initialize the ExtractorBase class.

        Sets up the configuration and logger. Initializes the database session if required.
        """
        self.conf = conf
        self.logger = setup_logger(self.__class__.__name__)
        self.logger.info(f"Initializing {self.__class__.__name__}")

        self.redis_conn = Redis(host='localhost', port=6379, db=0)
        self.process_queue = Queue(connection=self.redis_conn)
        self.data_queue = multiprocessing.Queue()

        self.queues = {}

        if session_required:
            self.session_init()

    def session_init(self):
        """
        Initialize the database session using DatabaseManager.

        Sets up the database connection and session using the DatabaseManager class.
        """
        self.logger.info("Initializing database session using DatabaseManager")
        db_manager = DatabaseManager(self.conf)
        self.engine = db_manager.get_engine()
        self.session = db_manager.get_session()

    @abstractmethod
    def set_tasks(self):
        pass

    def start(self):
        """
        Start the data extraction process.

        This abstract method should be implemented by all subclasses to define
        the specific data extraction logic for each bioinformatics data source.
        """
        self.set_tasks()
        # self.set_targets()
        # self.start_db_process()
        # self.fetch()
        # self.data_queue.put(None)  # Señal de terminación para el proceso de base de datos
        # self.db_process.join()  # Esperar a que el proceso de base de datos finalice

    def start_db_process(self):
        """Inicia el proceso que maneja la inserción en la base de datos."""
        self.db_process = multiprocessing.Process(target=self.add_to_db)
        self.db_process.start()

    def add_to_db(self):
        """
        Procesa elementos de la cola y los añade a la base de datos.
        """
        while True:
            data = self.data_queue.get()
            if data is None:  # Verificar si es la señal de terminación
                break
            try:
                self.store_entry(data)
            except Exception as e:
                self.logger.error(f"Error processing data: {str(e)}")
        self.session.close()
        self.logger.info("Database session closed and process ending.")

    # @abstractmethod
    # def set_targets(self):
    #     pass
    #
    # @abstractmethod
    # def fetch(self):
    #     pass
    #
    # @abstractmethod
    # def store_entry(self, data):
    #     pass
