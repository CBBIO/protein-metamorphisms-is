from abc import ABC


from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager


class BaseTaskInitializer(ABC):
    def __init__(self, conf, session_required=True):
        self.logger = setup_logger(self.__class__.__name__)
        self.logger.info(f"Initializing {self.__class__.__name__}")
        self.conf = conf

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

