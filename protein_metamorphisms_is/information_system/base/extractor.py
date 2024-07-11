import multiprocessing
import pickle
import threading
import time
import traceback

import pika
from pika import PlainCredentials
from abc import ABC, abstractmethod

from protein_metamorphisms_is.base.task import BaseTaskInitializer
from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.sql.base.database_manager import DatabaseManager

class ExtractorBase(BaseTaskInitializer, ABC):
    def __init__(self, conf, session_required=False):
        super().__init__(conf, session_required=True)
        self.stop_event = multiprocessing.Event()



