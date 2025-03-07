import importlib
import traceback

from protein_metamorphisms_is.sql.model.entities.embedding.sequence_embedding import (
    SequenceEmbeddingType,
    SequenceEmbedding,
)

from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.tasks.gpu import GPUTaskInitializer


class SequenceEmbeddingManager(GPUTaskInitializer):
    """
    Manages the sequence embedding process, including model loading, task enqueuing, and result storing.

    This class initializes GPU tasks, retrieves model configuration, and processes batches of sequences
    for embedding generation.

    Attributes:
        reference_attribute (str): Name of the attribute used as the reference for embedding (default: 'sequence').
        model_instances (dict): Dictionary of loaded models keyed by embedding type ID.
        tokenizer_instances (dict): Dictionary of loaded tokenizers keyed by embedding type ID.
        base_module_path (str): Base module path for dynamic imports of embedding tasks.
        batch_size (int): Number of sequences processed per batch. Defaults to 40.
        types (dict): Configuration dictionary for embedding types.
    """

    def __init__(self, conf):
        """
        Initializes the SequenceEmbeddingManager.

        :param conf: Configuration dictionary containing embedding parameters.
        :type conf: dict

        Example:
            >>> conf = {
            >>>     "embedding": {"batch_size": 50, "types": [1, 2]},
            >>>     "limit_execution": 100
            >>> }
            >>> manager = SequenceEmbeddingManager(conf)
        """
        super().__init__(conf)
        self.reference_attribute = 'sequence'
        self.model_instances = {}
        self.tokenizer_instances = {}
        self.base_module_path = 'protein_metamorphisms_is.operation.embedding.proccess.sequence'
        self.fetch_models_info()
        self.batch_size = self.conf['embedding'].get('batch_size', 40)

    def fetch_models_info(self):
        """
        Retrieves and initializes embedding models based on the database configuration.

        :raises sqlalchemy.exc.SQLAlchemyError: If there's an error querying the database.
        """
        self.session_init()
        embedding_types = self.session.query(SequenceEmbeddingType).all()
        self.session.close()
        del self.engine
        self.types = {}

        for type_obj in embedding_types:
            if type_obj.id in self.conf['embedding']['types']:
                module_name = f"{self.base_module_path}.{type_obj.task_name}"
                module = importlib.import_module(module_name)
                self.types[type_obj.id] = {
                    'module': module,
                    'model_name': type_obj.model_name,
                    'id': type_obj.id,
                    'task_name': type_obj.task_name,
                }

    def enqueue(self):
        """
        Enqueues tasks for sequence embedding processing.

        Splits all sequences into batches, checks for existing embeddings in the database,
        and prepares tasks for missing embeddings.

        :raises Exception: If an error occurs during the enqueue process.

        Example:
            >>> manager.enqueue()
        """
        try:
            self.logger.info("Starting embedding enqueue process.")
            self.session_init()
            sequences = self.session.query(Sequence).all()

            if self.conf['limit_execution']:
                sequences = sequences[:self.conf['limit_execution']]

            sequence_batches = [
                sequences[i: i + self.batch_size]
                for i in range(0, len(sequences), self.batch_size)
            ]

            for batch in sequence_batches:
                model_batches = {}
                for sequence in batch:
                    for type in self.types.values():
                        existing_embedding = self.session.query(SequenceEmbedding).filter_by(
                            sequence_id=sequence.id, embedding_type_id=type['id']
                        ).first()

                        if not existing_embedding:
                            task_data = {
                                'sequence': sequence.sequence,
                                'sequence_id': sequence.id,
                                'model_name': type['model_name'],
                                'embedding_type_id': type['id'],
                            }
                            model_batches.setdefault(type['id'], []).append(task_data)

                for model_type, batch_data in model_batches.items():
                    if batch_data:
                        self.publish_task(batch_data, model_type)
                        self.logger.info(
                            f"Published batch with {len(batch_data)} sequences to model type {model_type}."
                        )
            self.session.close()

        except Exception as e:
            self.logger.error(f"Error during enqueue process: {e}")
            raise

    def process(self, batch_data):
        """
        Processes a batch of sequences to generate embeddings.

        :param batch_data: List of dictionaries, each containing sequence data.
        :type batch_data: list[dict]
        :return: List of dictionaries with embedding results.
        :rtype: list[dict]
        :raises Exception: If there's an error during embedding generation.

        Example:
            >>> batch_data = [{"sequence": "ATCG", "sequence_id": 1, "embedding_type_id": 2}]
            >>> results = manager.process(batch_data)
        """
        try:
            embedding_type_id = batch_data[0]['embedding_type_id']
            model = self.model_instances[embedding_type_id]
            tokenizer = self.tokenizer_instances[embedding_type_id]
            module = self.types[embedding_type_id]['module']
            device = self.conf['embedding'].get('device', "cuda")

            embedding_records = module.embedding_task(
                batch_data, model, tokenizer, embedding_type_id=embedding_type_id, device=device
            )

            return embedding_records

        except Exception as e:
            self.logger.error(f"Error during embedding process: {e}\n{traceback.format_exc()}")
            raise

    def store_entry(self, records):
        """
        Stores embedding results in the database.

        :param records: List of embedding result dictionaries.
        :type records: list[dict]
        :raises RuntimeError: If an error occurs during database storage.

        Example:
            >>> records = [
            >>>     {"sequence_id": 1, "embedding_type_id": 2, "embedding": [0.1, 0.2], "shape": [2]}
            >>> ]
            >>> manager.store_entry(records)
        """
        session = self.session
        try:
            for record in records:
                embedding_entry = SequenceEmbedding(
                    sequence_id=record['sequence_id'],
                    embedding_type_id=record['embedding_type_id'],
                    embedding=record['embedding'],
                    shape=record['shape'],
                )
                session.add(embedding_entry)
            session.commit()

        except Exception as e:
            session.rollback()
            self.logger.error(f"Error during database storage: {e}")
            raise RuntimeError(f"Error storing entry: {e}")
