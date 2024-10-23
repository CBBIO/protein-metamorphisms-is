import importlib
from itertools import combinations
from operator import and_, or_

from sqlalchemy import func
from sqlalchemy.exc import SQLAlchemyError

from protein_metamorphisms_is.sql.model.entities.embedding.structure_3di import Structure3Di
from protein_metamorphisms_is.sql.model.entities.structure.state import State
from protein_metamorphisms_is.sql.model.entities.structure.structure import Structure
from protein_metamorphisms_is.sql.model.operational.clustering.cluster import SubclusterEntry, Subcluster
from protein_metamorphisms_is.sql.model.operational.structural_alignment.group import AlignmentGroup, \
    AlignmentGroupEntry
from protein_metamorphisms_is.sql.model.operational.structural_alignment.structural_alignment_result import \
    AlignmentResult
from protein_metamorphisms_is.sql.model.operational.structural_alignment.structural_alignment_type import \
    StructuralAlignmentType

from protein_metamorphisms_is.tasks.queue import QueueTaskInitializer
import logging


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
        self.types = {}
        base_module_path = 'protein_metamorphisms_is.operation.structural_alignment.tasks'

        for type_obj in structural_alignment_types:
            if type_obj.id in self.conf['structural_alignment']['types']:
                module_name = f"{base_module_path}.{type_obj.task_name}"
                module = importlib.import_module(module_name)
                self.types[type_obj.id] = module

    from sqlalchemy import or_, and_

    def enqueue(self):
        """
        Enqueue tasks for all pairs of representational subclusters within each cluster.
        The tasks are created based on the activated alignment types.
        """
        subclusters_query = self.session.query(
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
        ).join(
            Structure, Structure.id == State.structure_id
        ).filter(
            SubclusterEntry.is_representative == True
        ).all()

        # Group the subclusters by cluster_id
        clusters_dict = {}
        for entry in subclusters_query:
            clusters_dict.setdefault(entry.cluster_id, []).append(entry)

        # Process each cluster to generate pair combinations
        for cluster_id, subcluster_entries in clusters_dict.items():
            if len(subcluster_entries) < 2:
                continue
            pair_combinations = list(combinations(subcluster_entries, 2))

            # Enqueue tasks for each pair combination
            for subcluster_1, subcluster_2 in pair_combinations:
                # Enqueue tasks for each activated alignment type
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
            if result:
                self.store_entry(result)
            else:
                self.logger.error(result)
        except Exception as e:
            self.logger.error(f"Error processing task for alignment group {data.get('alignment_group_id')}: {str(e)}")

    def store_entry(self, record):
        """
        Store the results of the alignment in the AlignmentResult table and manage AlignmentGroup and AlignmentGroupEntry creation.
        """
        try:
            # Retrieve the necessary IDs from the task data
            cluster_id = record.get('cluster_id')
            subcluster_entry_1_id = record.get('subcluster_entry_1_id')
            subcluster_entry_2_id = record.get('subcluster_entry_2_id')

            # Check if an AlignmentGroup already exists for these entries
            existing_group = self.session.query(AlignmentGroup).join(AlignmentGroupEntry).filter(
                AlignmentGroupEntry.cluster_entry_id.in_([subcluster_entry_1_id, subcluster_entry_2_id])
            ).group_by(AlignmentGroup.id).having(func.count(AlignmentGroupEntry.id) == 2).first()

            if existing_group:
                alignment_group_id = existing_group.id
                self.logger.info(
                    f"Using existing alignment group {alignment_group_id} for entries {subcluster_entry_1_id} and {subcluster_entry_2_id}.")
            else:
                # Create a new AlignmentGroup
                alignment_group = AlignmentGroup()
                self.session.add(alignment_group)
                self.session.flush()  # Ensure the new AlignmentGroup is available for use
                alignment_group_id = alignment_group.id
                self.logger.info(f"Created a new alignment group with ID {alignment_group_id}.")

                # Create AlignmentGroupEntries for the pair
                entries = [
                    AlignmentGroupEntry(alignment_group_id=alignment_group_id, cluster_entry_id=subcluster_entry_1_id),
                    AlignmentGroupEntry(alignment_group_id=alignment_group_id, cluster_entry_id=subcluster_entry_2_id)
                ]
                self.session.add_all(entries)
                self.session.flush()  # Ensure the new AlignmentGroupEntries are staged
            print(alignment_group_id)
            # Check if the alignment result already exists
            existing_result = self.session.query(AlignmentResult).filter_by(
                alignment_group_id=alignment_group_id).first()

            if existing_result:
                # Update the existing result
                existing_result.ce_rms = record.get('ce_rms', existing_result.ce_rms)
                existing_result.tm_rms = record.get('tm_rms', existing_result.tm_rms)
                existing_result.tm_seq_id = record.get('tm_seq_id', existing_result.tm_seq_id)
                existing_result.tm_score_chain_1 = record.get('tm_score_chain_1', existing_result.tm_score_chain_1)
                existing_result.tm_score_chain_2 = record.get('tm_score_chain_2', existing_result.tm_score_chain_2)
                existing_result.fc_rms = record.get('fc_rms', existing_result.fc_rms)
                existing_result.fc_identity = record.get('fc_identity', existing_result.fc_identity)
                existing_result.fc_similarity = record.get('fc_similarity', existing_result.fc_similarity)
                existing_result.fc_score = record.get('fc_score', existing_result.fc_score)
                existing_result.fc_align_len = record.get('fc_align_len', existing_result.fc_align_len)
                self.logger.info(f"Updated existing alignment result for alignment group {alignment_group_id}.")
            else:
                # Insert a new alignment result if none exists
                new_result = AlignmentResult(
                    alignment_group_id=alignment_group_id,
                    ce_rms=record.get('ce_rms'),
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
                self.logger.info(f"Created a new alignment result for alignment group {alignment_group_id}.")

            # Commit the transaction
            self.session.commit()
            self.logger.info(f"Alignment results stored successfully for alignment group {alignment_group_id}.")

        except SQLAlchemyError as e:
            self.session.rollback()
            self.logger.error(f"Error storing alignment results: {str(e)}")





