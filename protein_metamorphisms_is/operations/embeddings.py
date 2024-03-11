import importlib
import multiprocessing
from datetime import datetime, timedelta
from multiprocessing import Pool

from sqlalchemy.orm import aliased


from protein_metamorphisms_is.operations.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import PDBChains, Cluster, PDBReference, StructuralAlignmentQueue, \
    StructuralAlignmentType, StructuralAlignmentResults, EmbeddingType


class EmbeddingManager(OperatorBase):
    """
    lorem ipsum

    Attributes:
        conf (dict): Configuration of the instance, including database connections and operational settings.
    """

    def __init__(self, conf):
        """
        Initializes an instance of `StructuralAlignmentManager` with configuration settings.

        Args:
            conf (dict): Configuration parameters, including database connections and operational settings.
        """
        super().__init__(conf)
        self.logger.info("Secuence Embedding Manager instance created.")

    def fetch_models_info(self):
        """
        Fetches and prepares alignment task modules based on the configuration.

        This method dynamically imports alignment task modules specified in the configuration and stores
        references to these modules in a dictionary for later use in the alignment process.
        """
        embedding_types = self.session.query(EmbeddingType).all()
        self.types = {}
        base_module_path = 'protein_metamorphisms_is.operations.embedding_tasks'

        for type_obj in embedding_types:
            if type_obj.id in self.conf['embedding']['types']:
                # Construye el nombre completo del módulo
                module_name = f"{base_module_path}.{type_obj.task_name}"
                # Importa dinámicamente el módulo usando importlib
                module = importlib.import_module(module_name)
                # Almacena la referencia al módulo en el diccionario self.types
                self.types[type_obj.id] = {'module': module, 'model_name' : type_obj.model_name}

        print(self.types)

    def start(self):
        """
        Begin the structural alignment process.

        This method manages the workflow of the alignment process, including loading clusters, executing alignments,
        and handling any exceptions encountered during the process. Progress and errors are logged appropriately.
        """
        try:
            self.logger.info("Starting structural alignment process.")
            chains = self.load_chains()
            self.fetch_models_info()

            for type in self.types.values():
                print(type)
                module, model = type['module'], type['model_name']
                module.embedding_task(self.session,chains,module,model)



        except Exception as e:
            self.logger.error(f"Error during structural alignment process: {e}")
            raise




