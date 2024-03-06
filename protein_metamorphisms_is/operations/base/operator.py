from abc import abstractmethod, ABC

import yaml

from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager
from protein_metamorphisms_is.sql.constants import handle_structural_complexity_levels, \
    handle_structural_alignment_types


class OperatorBase(ABC):
    def __init__(self, conf):
        self.conf = conf
        self.logger = setup_logger(self.__class__.__name__)
        self.logger.info(f"Initializing {self.__class__.__name__}")

        db_manager = DatabaseManager(conf)
        self.session = db_manager.get_session()
        constants = yaml.safe_load(open(conf['constants']))
        handle_structural_complexity_levels(self.session, constants)
        handle_structural_alignment_types(self.session, constants)

    @abstractmethod
    def start(self):
        """
        Start the operation over data process.

        This abstract method should be implemented by all subclasses to define
        the specific data operation logic for each bioinformatics data source.
        """
        pass
