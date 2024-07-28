from abc import ABC, abstractmethod

import yaml

from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager
from protein_metamorphisms_is.sql.constants import handle_structural_complexity_levels, \
    handle_structural_alignment_types, handle_prediction_methods, \
    handle_sequence_embedding_types, handle_structure_embedding_types


class BaseTaskInitializer(ABC):
    def __init__(self, conf, session_required=True):
        self.logger = setup_logger(self.__class__.__name__)
        self.logger.info(f"Initializing {self.__class__.__name__}")
        self.conf = conf



        if session_required:
            self.session_init()

            constants = yaml.safe_load(open(conf['constants']))
            handle_structural_complexity_levels(self.session, constants)
            handle_structural_alignment_types(self.session, constants)
            handle_sequence_embedding_types(self.session, constants)
            handle_structure_embedding_types(self.session, constants)
            handle_prediction_methods(self.session, constants)


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
        Start the operation over data process.

        This abstract method should be implemented by all subclasses to define
        the specific data operation logic for each bioinformatics data source.
        """
        pass

    @abstractmethod
    def process(self, target):
        """
        Process the given target.

        This abstract method should be implemented by all subclasses to define
        the specific processing logic for each bioinformatics data source.
        """
        pass

    @abstractmethod
    def store_entry(self, record):
        """
        Store the processed entry.

        This abstract method should be implemented by all subclasses to define
        the specific storage logic for processed data.
        """
        pass
