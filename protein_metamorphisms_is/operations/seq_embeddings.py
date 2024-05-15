import importlib
import multiprocessing
from datetime import datetime, timedelta
from multiprocessing import Pool

from sqlalchemy.orm import aliased

from protein_metamorphisms_is.operations.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import PDBChains, Cluster, PDBReference, StructuralAlignmentQueue, \
    StructuralAlignmentType, StructuralAlignmentResults, EmbeddingType, Sequence


class SequenceEmbeddingManager(OperatorBase):
    def __init__(self, conf):
        super().__init__(conf)
        self.logger.info("Sequence Embedding Manager instance created.")

    def fetch_models_info(self):
        embedding_types = self.session.query(EmbeddingType).all()
        self.types = {}
        base_module_path = 'protein_metamorphisms_is.operations.embedding_tasks'

        for type_obj in embedding_types:
            if type_obj.id in self.conf['embedding']['types']:
                module_name = f"{base_module_path}.{type_obj.task_name}"
                module = importlib.import_module(module_name)
                self.types[type_obj.id] = {'module': module, 'model_name': type_obj.model_name, 'id': type_obj.id}

    def start(self):
        try:
            self.logger.info("Starting embedding process.")
            sequences = self.session.query(Sequence).all()
            self.fetch_models_info()

            for type in self.types.values():
                module, model_name, embedding_type_id = type['module'], type['model_name'], type['id']
                module.embedding_task(self.session, sequences, model_name, embedding_type_id)

        except Exception as e:
            self.logger.error(f"Error during embedding process: {e}")
            raise
