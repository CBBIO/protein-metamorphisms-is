import multiprocessing
import queue
import time
from abc import ABC, abstractmethod

from redis import Redis
from rq import Queue, Worker

from protein_metamorphisms_is.base.task import BaseTaskInitializer
from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager


class ExtractorBase(BaseTaskInitializer):
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
        super().__init__(conf, session_required=True)
        self.redis_conn = Redis(host='localhost', port=6379, db=0)  # Asegúrate de configurar según tu entorno de Redis
        self.process_queue = Queue('uniprot_extractor',connection=self.redis_conn)
        self.data_queue = multiprocessing.Queue()

    def start_worker(self):
        print('h')
        self.logger.info("Starting RQ worker for UniProt data extraction")

        # Inicializar el trabajador en un proceso separado
        def worker_function():
            worker = Worker([self.process_queue], connection=self.redis_conn)
            worker.work()

        self.worker_process = multiprocessing.Process(target=worker_function)
        self.worker_process.start()

    def start(self):
        """
        Start the data extraction process.

        This abstract method should be implemented by all subclasses to define
        the specific data extraction logic for each bioinformatics data source.
        """
        self.start_worker()
        self.start_db_process()
        self.queue_in()
        self.fetch()
        self.data_queue.put(None)  # Señal de terminación para el proceso de base de datos
        self.db_process.join()  # Esperar a que el proceso de base de datos finalice

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
    # def queue_in(self):
    #     pass
    #
    # @abstractmethod
    # def extract(self):
    #     pass
    #
    # @abstractmethod
    # def insert_db(self, data):
    #     pass


    def queue_in(self):
        pass
