import importlib
import traceback

from protein_metamorphisms_is.base.gpu import GPUTaskInitializer
from protein_metamorphisms_is.sql.model import StructureEmbedding, Model, StructureEmbeddingType


class StructureEmbeddingManager(GPUTaskInitializer):
    def __init__(self, conf):
        super().__init__(conf)
        self.reference_attribute = 'model'
        self.model_instances = {}
        self.tokenizer_instances = {}
        self.fetch_models_info()
        self.batch_size = self.conf['embedding'].get('batch_size', 40)

    def fetch_models_info(self):
        self.session_init()
        embedding_types = self.session.query(StructureEmbeddingType).all()
        self.session.close()
        del self.engine
        self.types = {}
        base_module_path = 'protein_metamorphisms_is.operations.structure_embedding_tasks'

        for type_obj in embedding_types:
            if type_obj.id in self.conf['embedding']['types']:
                module_name = f"{base_module_path}.{type_obj.task_name}"
                module = importlib.import_module(module_name)
                self.types[type_obj.id] = {
                    'module': module,
                    'model_name': type_obj.model_name,
                    'id': type_obj.id,
                    'task_name': type_obj.task_name
                }

                model = module.load_model(type_obj.model_name)
                tokenizer = module.load_tokenizer(type_obj.model_name)
                self.model_instances[type_obj.id] = model
                self.tokenizer_instances[type_obj.id] = tokenizer

    def enqueue(self):
        try:
            self.logger.info("Starting embedding enqueue process.")
            self.session_init()
            models = self.session.query(Model).all()
            model_batches = [models[i:i + self.batch_size] for i in range(0, len(models), self.batch_size)]

            for batch in model_batches:
                model_tasks = {}
                for model in batch:
                    for type in self.types.values():
                        task_data = {
                            'model': model,
                            'model_id': model.id,
                            'model_name': type['model_name'],
                            'embedding_type_id': type['id']
                        }
                        if type['id'] not in model_tasks:
                            model_tasks[type['id']] = []
                        model_tasks[type['id']].append(task_data)

                for model_type, task_data in model_tasks.items():
                    self.publish_task(task_data, model_type)
                    self.logger.info(f"Published batch with {len(task_data)} models to model type {model_type}.")

            self.session.close()

        except Exception as e:
            self.logger.error(f"Error during enqueue process: {e}")
            raise

    def process(self, task_data):
        print(task_data)
        try:
            results = []

            for data in task_data:
                embedding_type_id = data['model']
                print(embedding_type_id)
                model = self.model_instances[data['embedding_type_id']]
                print('model',model)
                tokenizer = self.tokenizer_instances[data['embedding_type_id']]
                module = self.types[data['embedding_type_id']]['module']

                model_id = data['model']
                model_instance = self.session.query(Model).filter_by(id=data['model_id']).one()

                parser = MMCIFParser()
                structure = parser.get_structure(model_instance.model_id, file_path)

                embedding_records = module.embedding_task([structure], model, tokenizer)

                for record in embedding_records:
                    record['model_id'] = data['model_id']
                    record['embedding_type_id'] = embedding_type_id
                    results.append(record)
            return results
        except Exception as e:
            self.logger.error(f"Error during embedding process: {e}\n{traceback.format_exc()}")
            raise

    def store_entry(self, records):
        session = self.session
        try:
            for record in records:
                embedding_entry = StructureEmbedding(
                    model_id=record['model_id'],
                    embedding_type_id=record['embedding_type_id'],
                    embedding=record['embedding'],
                    shape=record['shape']
                )
                session.add(embedding_entry)
            session.commit()
        except Exception as e:
            session.rollback()
            raise RuntimeError(f"Error storing entry: {e}")
