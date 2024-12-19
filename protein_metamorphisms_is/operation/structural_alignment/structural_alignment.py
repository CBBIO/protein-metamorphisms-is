import importlib
from itertools import combinations
from sqlalchemy import func
from sqlalchemy.exc import SQLAlchemyError

from protein_metamorphisms_is.sql.model.entities.embedding.structure_3di import Structure3Di
from protein_metamorphisms_is.sql.model.entities.structure.state import State
from protein_metamorphisms_is.sql.model.operational.clustering.cluster import SubclusterEntry, Subcluster
from protein_metamorphisms_is.sql.model.operational.structural_alignment.group import AlignmentGroup, \
    AlignmentGroupEntry
from protein_metamorphisms_is.sql.model.operational.structural_alignment.result import \
    AlignmentResult
from protein_metamorphisms_is.sql.model.operational.structural_alignment.structural_alignment_type import \
    StructuralAlignmentType

from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer


class StructuralAlignmentManager(QueueTaskInitializer):
    """
    Manages the structural alignment process of representational subclusters using various alignment algorithms.
    """

    def __init__(self, conf):
        super().__init__(conf)
        self.logger.info("Structural Alignment Manager initialized.")
        self.fetch_tasks_info()

    def fetch_tasks_info(self):
        """
        Load alignment task modules dynamically based on configuration.
        """
        structural_alignment_types = self.session.query(StructuralAlignmentType).all()
        base_module_path = 'protein_metamorphisms_is.operation.structural_alignment.tasks'

        self.types = {
            type_obj.id: importlib.import_module(f"{base_module_path}.{type_obj.task_name}")
            for type_obj in structural_alignment_types
            if type_obj.id in self.conf['structural_alignment']['types']
        }

    def enqueue(self):
        """
        Enqueue tasks for all pairs of representational subclusters within each cluster.
        Only queries the database without modifications.
        """
        subclusters = self.session.query(
            SubclusterEntry.id.label("subcluster_entry_id"),
            State.file_path.label("file_path"),
            Subcluster.id.label("subcluster_id"),
            Subcluster.cluster_id.label("cluster_id")
        ).join(
            Subcluster, Subcluster.id == SubclusterEntry.subcluster_id
        ).join(
            Structure3Di, Structure3Di.id == SubclusterEntry.structure_3di_id
        ).join(
            State, State.id == Structure3Di.state_id
        ).filter(
            SubclusterEntry.is_representative.is_(True)
        ).all()

        clusters_dict = {}
        for entry in subclusters:
            clusters_dict.setdefault(entry.cluster_id, []).append(entry)

        # Process each cluster to generate pair combinations
        for cluster_id, subcluster_entries in clusters_dict.items():
            if len(subcluster_entries) < 2:
                continue
            self._enqueue_tasks_for_cluster(cluster_id, subcluster_entries)

        self.logger.info("All tasks enqueued without modifying the database.")

    def _enqueue_tasks_for_cluster(self, cluster_id, subcluster_entries):
        """
        Enqueue tasks for a specific cluster.
        """
        for subcluster_1, subcluster_2 in combinations(subcluster_entries, 2):
            if not self._check_if_pair_exists(subcluster_1.subcluster_entry_id, subcluster_2.subcluster_entry_id):
                self._enqueue_task_for_alignment_types(cluster_id, subcluster_1, subcluster_2)

    def _check_if_pair_exists(self, subcluster_entry_1_id, subcluster_entry_2_id):
        """
        Check if an alignment group already exists for the given pair of subclusters.
        Also checks for null values in the associated alignment result.
        """
        existing_group = self.session.query(AlignmentGroup).join(AlignmentGroupEntry).filter(
            AlignmentGroupEntry.subcluster_entry_id.in_([subcluster_entry_1_id, subcluster_entry_2_id])
        ).group_by(AlignmentGroup.id).having(func.count(AlignmentGroupEntry.id) == 2).first()

        if existing_group:
            # Check if the corresponding AlignmentResult has null values
            alignment_result = self.session.query(AlignmentResult).filter_by(
                alignment_group_id=existing_group.id
            ).first()

            # Ensure alignment_result is not None before accessing its attributes
            if alignment_result is None or any(
                    value is None for value in [
                        alignment_result.ce_rms,
                        alignment_result.tm_rms,
                        alignment_result.tm_seq_id,
                        alignment_result.tm_score_chain_1,
                        alignment_result.tm_score_chain_2,
                        alignment_result.fc_rms,
                        alignment_result.fc_identity,
                        alignment_result.fc_similarity,
                        alignment_result.fc_score,
                        alignment_result.fc_align_len
                    ]
            ):
                self.logger.info(f"Alignment group {existing_group.id} exists but has null values. Allowing enqueue.")
                return False  # Allow enqueuing if values are null or result doesn't exist

            self.logger.info(
                f"Alignment group already exists for entries {subcluster_entry_1_id} and {subcluster_entry_2_id}. Skipping enqueue.")
            return True  # Skip enqueuing if group exists without nulls

        return False  # Return False if no existing group

    def _enqueue_task_for_alignment_types(self, cluster_id, subcluster_1, subcluster_2):
        """
        Enqueue tasks for each activated alignment type.
        """
        for alignment_type_id in self.conf['structural_alignment']['types']:
            task_data = {
                'cluster_id': cluster_id,
                'subcluster_entry_1_id': subcluster_1.subcluster_entry_id,
                'subcluster_entry_2_id': subcluster_2.subcluster_entry_id,
                'subcluster_1_file_path': subcluster_1.file_path,
                'subcluster_2_file_path': subcluster_2.file_path,
                'alignment_type_id': alignment_type_id
            }
            self.publish_task(task_data)
            self.logger.info(
                f"Task enqueued for subclusters {subcluster_1.subcluster_entry_id} and {subcluster_2.subcluster_entry_id} in cluster {cluster_id}")

    def process(self, data):
        """
        Process each task, executing the appropriate alignment algorithm for the alignment group.
        """
        try:
            self.logger.info(f"Processing task: {data}")
            alignment_type_id = data.get('alignment_type_id')

            if not alignment_type_id:
                raise ValueError("Missing required fields in task data: 'alignment_type_id'")

            align_task_module = self.types.get(alignment_type_id)
            if not align_task_module:
                raise ValueError(f"Alignment module for type_id {alignment_type_id} not found")

            result = align_task_module.align_task(data, self.conf, self.logger)
            return result if result else self.logger.error("No result returned from alignment task.")
        except Exception as e:
            self.logger.error(f"Error processing task for alignment group {data.get('alignment_group_id')}: {str(e)}")

    def store_entry(self, record):
        """
        Store the results of the alignment in the AlignmentResult table and manage AlignmentGroup and AlignmentGroupEntry creation.
        """
        try:
            subcluster_entry_1_id = record.get('subcluster_entry_1_id')
            subcluster_entry_2_id = record.get('subcluster_entry_2_id')

            if not subcluster_entry_1_id or not subcluster_entry_2_id:
                self.logger.error(f"Missing subcluster entry IDs: {subcluster_entry_1_id}, {subcluster_entry_2_id}")
                return

            # Get or create an alignment group and add entries to it
            alignment_group_id = self._get_or_create_alignment_group(subcluster_entry_1_id, subcluster_entry_2_id)

            if alignment_group_id is None:
                return  # Skip processing if the group is not valid

            self._create_or_update_alignment_result(alignment_group_id, record)

        except SQLAlchemyError as e:
            self.session.rollback()
            self.logger.error(f"Error storing alignment results: {str(e)}")

    def _get_or_create_alignment_group(self, subcluster_entry_1_id, subcluster_entry_2_id):
        """
        Get or create an AlignmentGroup for the specified subcluster entries.
        """
        existing_group = self.session.query(AlignmentGroup).join(AlignmentGroupEntry).filter(
            AlignmentGroupEntry.subcluster_entry_id.in_([subcluster_entry_1_id, subcluster_entry_2_id])
        ).group_by(AlignmentGroup.id).having(func.count(AlignmentGroupEntry.id) == 2).first()

        if existing_group:
            return existing_group.id

        # Create a new AlignmentGroup
        alignment_group = AlignmentGroup()
        self.session.add(alignment_group)
        self.session.flush()  # Ensure the new AlignmentGroup is available for use
        self.logger.info(f"Created a new alignment group with ID {alignment_group.id}.")

        # Create entries for the new alignment group
        self._create_alignment_group_entries(alignment_group.id, [subcluster_entry_1_id, subcluster_entry_2_id])

        return alignment_group.id

    def _create_alignment_group_entries(self, alignment_group_id, subcluster_entry_ids):
        """
        Create AlignmentGroupEntry records for the specified alignment group.
        """
        for subcluster_entry_id in subcluster_entry_ids:
            new_entry = AlignmentGroupEntry(
                alignment_group_id=alignment_group_id,
                subcluster_entry_id=subcluster_entry_id
            )
            self.session.add(new_entry)
            self.logger.info(
                f"Added new entry for subcluster_entry_id {subcluster_entry_id} in alignment group {alignment_group_id}.")

    def _create_or_update_alignment_result(self, alignment_group_id, record):
        """
        Create or update the alignment result for the specified alignment group.
        """
        existing_result = self.session.query(AlignmentResult).filter_by(
            alignment_group_id=alignment_group_id
        ).first()

        alignment_result_data = {
            'alignment_group_id': alignment_group_id,
            'ce_rms': record.get('ce_rms'),
            'tm_rms': record.get('tm_rms'),
            'tm_seq_id': record.get('tm_seq_id'),
            'tm_score_chain_1': record.get('tm_score_chain_1'),
            'tm_score_chain_2': record.get('tm_score_chain_2'),
            'fc_rms': record.get('fc_rms'),
            'fc_identity': record.get('fc_identity'),
            'fc_similarity': record.get('fc_similarity'),
            'fc_score': record.get('fc_score'),
            'fc_align_len': record.get('fc_align_len')
        }

        if existing_result:
            for key, value in alignment_result_data.items():
                if value is not None:
                    setattr(existing_result, key, value)
            self.logger.info(f"Updated existing alignment result for alignment group {alignment_group_id}.")
        else:
            new_result = AlignmentResult(**alignment_result_data)
            self.session.add(new_result)
            self.logger.info(f"Created a new alignment result for alignment group {alignment_group_id}.")

        self.session.commit()
        self.logger.info(f"Alignment results stored successfully for alignment group {alignment_group_id}.")
