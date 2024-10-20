import importlib
from itertools import combinations
from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
from protein_metamorphisms_is.sql.model.model import SubclusterEntry, Subcluster, StructureEmbedding, Model, Structure, \
    StructuralAlignmentResults, StructuralAlignmentType
import logging


class StructuralAlignmentManager(QueueTaskInitializer):
    """
    Manages the structural alignment process of representational subclusters using various alignment algorithms.
    """

    def __init__(self, conf):
        super().__init__(conf)
        self.logger = logging.getLogger(__name__)
        self.logger.info("Structural Alignment Manager initialized.")
        self.fetch_tasks_info()

    def fetch_tasks_info(self):
        """
        Load alignment task modules dynamically based on configuration.
        """
        structural_alignment_types = self.session.query(StructuralAlignmentType).all()
        self.types = {}
        base_module_path = 'protein_metamorphisms_is.operation.structural_alignment.tasks'

        for type_obj in structural_alignment_types:
            if type_obj.id in self.conf['structural_alignment']['types']:
                module_name = f"{base_module_path}.{type_obj.task_name}"
                module = importlib.import_module(module_name)
                self.types[type_obj.id] = module

    def enqueue(self):
        """
        Enqueue tasks for all pairs of representational subclusters within each cluster, including file paths.
        The task is created for each activated alignment type.
        """
        subclusters_query = self.session.query(
            SubclusterEntry.id.label("subcluster_entry_id"),
            Model.file_path.label("file_path"),
            Subcluster.id.label("subcluster_id"),
            Subcluster.cluster_id.label("cluster_id")
        ).join(
            Subcluster, Subcluster.id == SubclusterEntry.subcluster_id
        ).join(
            StructureEmbedding, StructureEmbedding.id == SubclusterEntry.structure_embedding_id
        ).join(
            Model, Model.id == StructureEmbedding.model_id
        ).join(
            Structure, Structure.id == Model.structure_id
        ).filter(
            SubclusterEntry.is_representative == True
        ).all()

        # Agrupamos los subclusters por cluster_id
        clusters_dict = {}
        for entry in subclusters_query:
            if entry.cluster_id not in clusters_dict:
                clusters_dict[entry.cluster_id] = []
            clusters_dict[entry.cluster_id].append(entry)

        # Publicación de tareas con las combinaciones de pares para cada tipo de alineamiento activado
        total = 0
        for cluster_id, subcluster_entries in clusters_dict.items():
            pair_combinations = list(combinations(subcluster_entries, 2))
            total+= len(pair_combinations)
            for subcluster_1, subcluster_2 in pair_combinations:
                for alignment_type_id in self.conf['structural_alignment']['types']:  # Iterar sobre los tipos activados
                    task_data = {
                        'cluster_id': cluster_id,
                        'subcluster_entry_1_id': subcluster_1.subcluster_entry_id,  # Acceso correcto al atributo
                        'subcluster_1_file_path': subcluster_1.file_path,  # Acceso correcto al atributo
                        'subcluster_entry_2_id': subcluster_2.subcluster_entry_id,  # Acceso correcto al atributo
                        'subcluster_2_file_path': subcluster_2.file_path,  # Acceso correcto al atributo
                        'subcluster_1_id': subcluster_1.subcluster_id,  # Añadido el subcluster_1_id
                        'subcluster_2_id': subcluster_2.subcluster_id,  # Añadido el subcluster_2_id
                        'alignment_type_id': alignment_type_id  # Publicar la tarea con el tipo de alineamiento
                    }
                    # Publicar tarea en RabbitMQ
                    self.publish_task(task_data)
                    self.logger.info(
                        f"Publishing task for cluster {cluster_id} with subclusters {subcluster_1.subcluster_id} and {subcluster_2.subcluster_id} for alignment type {alignment_type_id}")
        self.logger.info("All alignment tasks with file paths enqueued.")

    def process(self, data):
        """
        Process each task, executing the appropriate alignment algorithm for the subcluster pair.
        """
        try:
            self.logger.info(f"Processing task: {data}")
            subcluster_1_id = data['subcluster_entry_1_id']
            subcluster_2_id = data['subcluster_entry_2_id']
            alignment_type_id = data['alignment_type_id']

            # Fetch the appropriate alignment module
            align_task_module = self.types.get(alignment_type_id)

            if not align_task_module:
                raise ValueError(f"Alignment module for type_id {alignment_type_id} not found")

            # Execute the alignment task
            result = align_task_module.align_task(data, self.conf)

            self.logger.info(
                f"Task processed successfully for subclusters {subcluster_1_id} and {subcluster_2_id} using alignment type {alignment_type_id}")
            return result

        except Exception as e:
            # Logging the error in case of failure
            self.logger.error(
                f"Error processing task for subclusters {data.get('subcluster_entry_1_id')} and {data.get('subcluster_entry_2_id')}: {str(e)}")
            return None

    def store_entry(self, record):
        """
        Store the results of the alignment in the StructuralAlignmentResults table.
        If an entry with the same cluster_id, subcluster_1_id, and subcluster_2_id exists,
        it will update only the fields present in the record, leaving other fields unchanged.
        """
        try:
            # Buscar si ya existe una entrada con los mismos cluster_id, subcluster_1_id y subcluster_2_id
            existing_entry = self.session.query(StructuralAlignmentResults).filter_by(
                cluster_id=record['cluster_id'],
                subcluster_1_id=record['subcluster_entry_1_id'],
                subcluster_2_id=record['subcluster_entry_2_id']
            ).first()

            print(record)

            # Si existe, actualizamos solo los campos presentes en el record
            if existing_entry:
                if 'ce_rms' in record:
                    existing_entry.ce_rms = record['ce_rms']
                if 'tm_rms' in record:
                    existing_entry.tm_rms = record['tm_rms']
                if 'tm_seq_id' in record:
                    existing_entry.tm_seq_id = record['tm_seq_id']
                if 'tm_score_chain_1' in record:
                    existing_entry.tm_score_chain_1 = record['tm_score_chain_1']
                if 'tm_score_chain_2' in record:
                    existing_entry.tm_score_chain_2 = record['tm_score_chain_2']
                if 'fc_rms' in record:
                    existing_entry.fc_rms = record['fc_rms']
                if 'fc_identity' in record:
                    existing_entry.fc_identity = record['fc_identity']
                if 'fc_similarity' in record:
                    existing_entry.fc_similarity = record['fc_similarity']
                if 'fc_score' in record:
                    existing_entry.fc_score = record['fc_score']
                if 'fc_align_len' in record:
                    existing_entry.fc_align_len = record['fc_align_len']

                self.logger.info(
                    f"Updated alignment result for cluster {record['cluster_id']} with subclusters {record['subcluster_entry_1_id']} and {record['subcluster_entry_2_id']}")
            else:
                # Si no existe, creamos una nueva entrada
                new_result = StructuralAlignmentResults(
                    cluster_id=record['cluster_id'],
                    subcluster_1_id=record['subcluster_entry_1_id'],
                    subcluster_2_id=record['subcluster_entry_2_id'],
                    ce_rms=record.get('ce_rms'),  # Usamos get para evitar None si no está presente
                    tm_rms=record.get('tm_rms'),
                    tm_seq_id=record.get('tm_seq_id'),
                    tm_score_chain_1=record.get('tm_score_chain_1'),
                    tm_score_chain_2=record.get('tm_score_chain_2'),
                    fc_rms=record.get('fc_rms'),
                    fc_identity=record.get('fc_identity'),
                    fc_similarity=record.get('fc_similarity'),
                    fc_score=record.get('fc_score'),
                    fc_align_len=record.get('fc_align_len')
                )
                self.session.add(new_result)
                self.logger.info(
                    f"Stored new alignment result for cluster {record['cluster_id']} with subclusters {record['subcluster_entry_1_id']} and {record['subcluster_entry_2_id']}")

            # Confirmamos los cambios
            self.session.commit()

        except Exception as e:
            self.session.rollback()  # Deshacer cambios en caso de error
            self.logger.error(f"Error storing or updating alignment result: {str(e)}")
