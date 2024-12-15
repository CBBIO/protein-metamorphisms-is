import importlib
import os
import traceback

from Bio import SeqIO

from protein_metamorphisms_is.operation.embedding.sequence_embedding import SequenceEmbeddingManager
from protein_metamorphisms_is.sql.model.entities.embedding.sequence_embedding import SequenceEmbeddingType, \
    SequenceEmbedding
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.tasks.gpu import GPUTaskInitializer
import h5py


class SequenceEmbedder(SequenceEmbeddingManager):
    def __init__(self, conf):
        super().__init__(conf)
        self.reference_attribute = 'sequence_embedder_from_fasta'
        self.model_instances = {}
        self.tokenizer_instances = {}
        self.base_module_path = 'protein_metamorphisms_is.operation.embedding.proccess.sequence'
        self.fetch_models_info()
        self.batch_size = self.conf['embedding'].get('batch_size', 40)  # Add batch size configuration
        self.fasta_path = conf.get('fantasia_input_fasta')
        self.output_csv = conf.get("fantasia_output_csv")
        self.length_filter = conf.get('length_filter', None)
        self.output_h5 = conf.get("fantasia_output_h5", "embeddings_results.h5")

        self.results = []

    def enqueue(self):
        try:
            self.logger.info("Starting embedding enqueue process.")
            sequences = []

            if self.conf.get('redundancy_filter'):
                os.system(
                    f"cd-hit -i {self.fasta_path} -o {self.conf.get('redundancy_file')} -c {self.conf.get('redundancy_filter')}")
                print(f"cd-hit -i {self.fasta_path} -o {self.conf.get('redundancy_file')} -c {self.conf.get('redundancy_filter')}")

            for record in SeqIO.parse(self.fasta_path, "fasta"):
                if self.length_filter and len(record.seq) > self.length_filter:
                    continue
                sequences.append(record)

            sequence_batches = [sequences[i:i + self.batch_size] for i in range(0, len(sequences), self.batch_size)]

            for batch in sequence_batches:
                model_batches = {}
                for sequence in batch:
                    for type in self.types.values():
                        task_data = {
                            'sequence': str(sequence.seq),
                            'accession': sequence.id,  # Usar el ID del FASTA como accession
                            'model_name': type['model_name'],
                            'embedding_type_id': type['id']
                        }

                        if type['id'] not in model_batches:
                            model_batches[type['id']] = []
                        model_batches[type['id']].append(task_data)

                for model_type, batch_data in model_batches.items():
                    self.publish_task(batch_data, model_type)
                    self.logger.info(
                        f"Published batch with {len(batch_data)} sequences to model type {model_type}.")

        except Exception as e:
            self.logger.error(f"Error during enqueue process: {e}")
            raise

    def process(self, task_data):
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
                    record['embedding_type_id'] = embedding_type_id
                    record['accession'] = data['accession']  # Propagar el accession
                    results.append(record)
            return results
        except Exception as e:
            self.logger.error(f"Error during embedding process: {e}\n{traceback.format_exc()}")
            raise

    def store_entry(self, results):
        try:
            # Verificar y crear el directorio si no existe
            output_dir = os.path.dirname(self.output_h5)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                self.logger.info(f"Created directory: {output_dir}")

            # Abrir el archivo en modo append ("a")
            with h5py.File(self.output_h5, "a") as h5file:
                for record in results:
                    accession = record['accession']
                    embedding_type_id = record['embedding_type_id']
                    print(embedding_type_id)

                    # Crear un grupo para el accession si no existe
                    accession_group = h5file.require_group(f"accession_{accession}")

                    # Crear un subgrupo para el embedding_type_id si no existe
                    type_group = accession_group.require_group(f"type_{embedding_type_id}")

                    # Evitar sobrescribir si el dataset ya existe
                    if "embedding" in type_group:
                        self.logger.warning(
                            f"Embedding for type {embedding_type_id} already exists in accession {accession}. Skipping.")
                        continue

                    # Crear el dataset
                    type_group.create_dataset("embedding", data=record['embedding'])
                    type_group.attrs['shape'] = record['shape']
                    self.logger.info(f"Stored embedding for accession {accession}, type {embedding_type_id}.")
        except Exception as e:
            self.logger.error(f"Error storing results in HDF5: {e}")
            raise

