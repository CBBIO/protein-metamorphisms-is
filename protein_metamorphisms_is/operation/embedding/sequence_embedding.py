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
        reference_attribute (str): Reference attribute name used for embedding.
        model_instances (dict): A dictionary storing loaded models.
        tokenizer_instances (dict): A dictionary storing loaded tokenizers.
        base_module_path (str): Base module path for dynamic model imports.
        batch_size (int): Batch size for processing sequences, defaults to 40.
        types (dict): Dictionary containing model and task configurations.
    """

    def __init__(self, conf):
        """
        Initializes the SequenceEmbeddingManager.

        Args:
            conf (dict): Configuration dictionary containing embedding parameters.
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
        Retrieves and initializes embedding models based on configuration.

        Queries the `SequenceEmbeddingType` table to fetch available embedding models.
        Modules are dynamically imported and stored in the `types` attribute.
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

        Fetches all sequences from the database, batches them, and checks for existing embeddings.
        If no existing embeddings are found, tasks are prepared and published for processing.

        Raises:
            Exception: If an error occurs during the enqueue process.
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

    def process(self, task_data):
        """
        Processes a batch of task data for embedding generation.

        Args:
            task_data (list[dict]): List of task data dictionaries, each containing:
                - sequence (str): Input sequence.
                - sequence_id (int): Sequence identifier.
                - embedding_type_id (int): Embedding type identifier.

        Returns:
            list[dict]: List of records with embedding results.

        Raises:
            Exception: If an error occurs during embedding processing.
        """
        try:
            results = []

            for data in task_data:
                embedding_type_id = data['embedding_type_id']
                model = self.model_instances[embedding_type_id]
                tokenizer = self.tokenizer_instances[embedding_type_id]
                module = self.types[embedding_type_id]['module']

                sequence = data['sequence']
                embedding_records = module.embedding_task([sequence], model, tokenizer)

                for record in embedding_records:
                    record['sequence_id'] = data['sequence_id']
                    record['embedding_type_id'] = embedding_type_id
                    results.append(record)
            return results

        except Exception as e:
            self.logger.error(f"Error during embedding process: {e}\n{traceback.format_exc()}")
            raise

    def store_entry(self, records):
        """
        Stores embedding records in the database.

        Args:
            records (list[dict]): List of embedding records, each containing:
                - sequence_id (int): Sequence identifier.
                - embedding_type_id (int): Embedding type identifier.
                - embedding: The embedding result.
                - shape: Shape of the embedding.

        Raises:
            RuntimeError: If an error occurs during database storage.
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
            raise RuntimeError(f"Error storing entry: {e}")
