from abc import abstractmethod, ABC

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from protein_data_handler.helpers.logger.logger import setup_logger
from protein_data_handler.operations.base.operator import OperatorBase
from protein_data_handler.sql.model import Base


class StructuralAlignmentBase(OperatorBase):
    def __init__(self, conf, session_required):
        self.conf = conf
        self.logger = setup_logger(self.__class__.__name__)
        self.logger.info(f"Initializing {self.__class__.__name__}")

        if session_required:
            self.session_init()


