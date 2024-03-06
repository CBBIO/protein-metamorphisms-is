from abc import ABC, abstractmethod

from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager


class ExtractorBase(ABC):
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
    def start(self):
        """
        Start the data extraction process.

        This abstract method should be implemented by all subclasses to define
        the specific data extraction logic for each bioinformatics data source.
        """
        pass
