import importlib
import logging
import multiprocessing
import os
from datetime import datetime, timedelta
from multiprocessing import Pool

from Bio.PDB import CEAligner, MMCIFParser
from sqlalchemy.orm import aliased

from protein_data_handler.operations.structural_alignment_tasks import universal, combinatorial_extension
from protein_data_handler.operations.base.operator import OperatorBase
from protein_data_handler.sql.model import PDBChains, Cluster, PDBReference, StructuralAlignmentQueue, \
    StructuralAlignmentType


class StructuralAlignmentManager(OperatorBase):
    """
    """

    def __init__(self, conf):
        """
        Initializes the CEAlign object with configuration settings.

        Args:
            conf (dict): Configuration parameters, including database connections and operational settings.
        """
        super().__init__(conf)
        self.logger.info("CEAlign instance created with configuration.")

    def fetch_tasks_info(self):
        structural_alignment_types = self.session.query(StructuralAlignmentType).all()
        self.types = {}
        base_module_path = 'protein_data_handler.operations.structural_alignment_tasks'

        for type_obj in structural_alignment_types:
            if type_obj.id in self.conf['structural_alignment']['types']:
                # Construye el nombre completo del módulo
                module_name = f"{base_module_path}.{type_obj.task_name}"
                # Importa dinámicamente el módulo usando importlib
                module = importlib.import_module(module_name)
                # Almacena la referencia al módulo en el diccionario self.types
                self.types[type_obj.id] = module

        print(self.types)

    def start(self):
        """
        Begin the structural alignment process using the CE algorithm.

        This method manages the workflow of the alignment process, including loading clusters, executing alignments,
        and handling any exceptions encountered during the process. Progress and errors are logged appropriately.
        """
        try:
            self.logger.info("Starting structural alignment process.")
            self.fetch_tasks_info()
            self.get_update_queue()

            while not self.check_empty_queue():
                queue_items = self.fetch_queue_items()
                self.execute_aligns(queue_items)

        except Exception as e:
            self.logger.error(f"Error during structural alignment process: {e}")
            raise

    def get_update_queue(self):
        for type_entry in self.types.items():
            type_id = type_entry[0]
            type_task_name = type_entry[1]
            self.logger.info(f"Type: {type_task_name} - Loading cluster entries from database.")

            queued_cluster_ids = self.session.query(StructuralAlignmentQueue.cluster_entry_id) \
                .filter(StructuralAlignmentQueue.alignment_type_id == type_id) \
                .all()
            self.logger.info(f"Found a total of {len(queued_cluster_ids)} clusters already in queue.")
            queued_cluster_ids = [c[0] for c in queued_cluster_ids]

            clusters_not_queued = self.session.query(Cluster).filter(
                Cluster.id.notin_(queued_cluster_ids),
                Cluster.is_representative == False
            ).all()

            self.logger.info(f"Found {len(clusters_not_queued)} clusters not in queue, adding to queue.")

            for cluster in clusters_not_queued:
                new_queue_item = StructuralAlignmentQueue(
                    cluster_entry_id=cluster.id,
                    state=0,
                    alignment_type_id=type_id,
                    retry_count=0
                )
                self.session.add(new_queue_item)
            #
            self.logger.info("New clusters added to queue. Checking for stale entries.")

            retry_timeout_date = datetime.now() - timedelta(minutes=self.conf.get('retry_timeout', 60))

            stale_entries = self.session.query(StructuralAlignmentQueue).filter(
                StructuralAlignmentQueue.state.in_([1, 3]),
                StructuralAlignmentQueue.alignment_type_id == type_id,
                StructuralAlignmentQueue.updated_at < retry_timeout_date
            ).all()
            self.logger.info(f"Found {len(stale_entries)} stale queue entries, resetting state to 'pending'.")

            for entry in stale_entries:
                entry.state = 0
                entry.retry_count += 1
                self.session.add(entry)

            self.session.commit()
            self.logger.info(f"task: {type_task_name} Queue update process completed.")

    def check_empty_queue(self):
        """
        Checks whether the alignment queue is empty or not.

        This method evaluates the current state of the CE alignment queue (`CEAlignQueue`) to determine if there are
        any pending or error entries that are still awaiting processing. It counts the number of entries in a 'pending'
        (state 0) or 'error' (state 3) state, with a retry count less than the configured limit.

        Returns:
            bool: Returns True if the queue is empty (i.e., no pending or error entries awaiting processing).
                  Returns False if there are entries in the queue that need to be processed.
        """
        queue_count = self.session.query(StructuralAlignmentQueue).filter(
            StructuralAlignmentQueue.state.in_([0, 3]),
            StructuralAlignmentQueue.alignment_type_id.in_(self.types.keys()),
            StructuralAlignmentQueue.retry_count < self.conf.get('retry_count', 5)
        ).count()

        is_empty = queue_count == 0

        if is_empty:
            self.logger.info("The alignment queue is empty.")
        else:
            self.logger.info(f"The alignment queue has {queue_count} items.")

        return is_empty

    def fetch_queue_items(self):
        batch_size = self.conf['structural_alignment']['batch_size']

        # Subconsulta para encontrar la estructura representativa para cada cluster
        representative_structures_subquery = self.session.query(
            Cluster.cluster_id,
            PDBReference.pdb_id.label('rep_pdb_id'),  # Cambiado para obtener pdb_id directamente
            PDBChains.chains.label('rep_chains'),
            PDBChains.model.label('rep_model')
        ).join(
            PDBChains, Cluster.pdb_chain_id == PDBChains.id
        ).join(
            PDBReference, PDBChains.pdb_reference_id == PDBReference.id  # Join adicional para acceder a PDBReference
        ).filter(
            Cluster.is_representative == True
        ).subquery()

        # Alias para hacer más fácil la referencia a la subconsulta
        rep_structures = aliased(representative_structures_subquery, name='rep_structures')

        # Realiza la consulta principal con un join a la subconsulta
        queue_items = self.session.query(
            StructuralAlignmentQueue.id,
            StructuralAlignmentQueue.alignment_type_id,
            Cluster.id.label("cluster_id"),
            PDBReference.pdb_id,
            PDBChains.chains,
            PDBChains.model,
            rep_structures.c.rep_pdb_id,
            rep_structures.c.rep_chains,
            rep_structures.c.rep_model
        ).join(
            Cluster, StructuralAlignmentQueue.cluster_entry_id == Cluster.id
        ).join(
            PDBChains, Cluster.pdb_chain_id == PDBChains.id
        ).join(
            PDBReference, PDBChains.pdb_reference_id == PDBReference.id
        ).outerjoin(
            rep_structures, Cluster.cluster_id == rep_structures.c.cluster_id
        ).filter(
            StructuralAlignmentQueue.state.in_([0, 3]),
            StructuralAlignmentQueue.alignment_type_id.in_(self.types.keys()),
            StructuralAlignmentQueue.retry_count < self.conf.get('retry_count', 5)
        ).order_by(
            StructuralAlignmentQueue.updated_at.desc()
        ).limit(
            batch_size
        ).all()

        return queue_items

    def execute_aligns(self, queue_items):
        num_workers = self.conf.get('max_workers', 4)
        timeout = self.conf.get('task_timeout', 1)
        self.logger.info(f"Performing CE alignment with {num_workers} workers.")

        results = []

        with Pool(num_workers) as pool:
            print(queue_items[0]._asdict())
            tasks = [pool.apply_async(self.types[item.alignment_type_id].align_task, args=(item, self.conf)) for item
                     in queue_items]

            for task in tasks:
                try:
                    result = task.get(timeout=timeout)
                    results.append(result)
                except multiprocessing.TimeoutError:
                    self.logger.error("Timeout error: La tarea excedió el tiempo límite establecido")

        self.logger.info(f"CE alignment completed for {len(results)} entries.")

        # self.insert_results(results)
        self.logger.info("CE alignment results inserted into the database successfully.")

    def insert_results(self, results):
        """
        Inserts the results of CE structural alignment tasks into the database and updates the alignment queue.

        This method processes the results obtained from the CE alignment tasks. For each result, it either stores the
        alignment data in the database or handles any errors that occurred during the alignment process. The method also
        updates the state of each queue item based on the outcome of its corresponding alignment task.

        Args:
            results (list): A list of dictionaries, each containing either the alignment results or error information.
                            The dictionaries include the cluster entry ID and, depending on success or failure,
                            either the RMS (Root Mean Square deviation) value or an error message.

        Process:
            - For successful alignments, the method stores the RMS value in the `CEAlignResults` table.
            - For failed alignments, indicated by an 'error_message' key in the result, the method logs the error and
              updates the queue item's state to 'error' (state 3), incrementing the retry count.
            - The method then updates the queue item's state to 'completed' (state 2) for successful alignments.

        Note:
            - The state of each queue item in `CEAlignQueue` is updated to reflect the result of the alignment task.
            - The states are defined as 0 for 'pending', 1 for 'in process', 2 for 'completed', and 3 for 'error'.
        """
        self.logger.info("Processing CE alignment results.")

        for result in results:
            if 'error_message' in result.keys():
                self.logger.warn(f"Error in alignment task: {result}")

                queue_item = self.session.query(CEAlignQueue).filter_by(
                    cluster_entry_id=result['cluster_entry_id']).first()
                if queue_item:
                    queue_item.state = 3  # Error state
                    queue_item.retry_count += 1
                    queue_item.error_message = str(result['error_message'])
                    self.session.add(queue_item)
            else:
                new_result = CEAlignResults(
                    cluster_entry_id=result['cluster_entry_id'],
                    rms=result['rms']
                )
                self.session.add(new_result)

                queue_item = self.session.query(CEAlignQueue).filter_by(
                    cluster_entry_id=result['cluster_entry_id']).first()
                if queue_item:
                    queue_item.state = 2
                    self.session.add(queue_item)

        self.session.commit()
        self.logger.info("CE alignment results processed and queue updated.")


#
# def align_task(alignment_entry):
#     """
#     Performs the alignment task for a given pair of protein structures.
#
#     This function aligns a target protein structure with a representative structure using the CE algorithm. It logs the
#     process and handles any exceptions that occur during the alignment. The function is designed to be run asynchronously
#     in a multiprocessing environment.
#
#     Args:
#         alignment_entry (dict): A dictionary containing the representative and target protein structures' information and paths.
#
#     Returns:
#         dict: A dictionary with the alignment result or error message, if any.
#     """
#     align_task_logger = logging.getLogger("align_task")
#     try:
#         align_task_logger.info("Aligning structures...")
#         representative = alignment_entry['representative']
#         target = alignment_entry['target']
#         pdb_chains_path = alignment_entry['pdb_chains_path']
#         parser = MMCIFParser()
#         representative_name = f"{representative['pdb_id']}_{representative['chains']}_{representative['model']}"
#         representative_structure_path = os.path.join(pdb_chains_path, f"{representative_name}.cif")
#         representative_structure = parser.get_structure(representative_name, representative_structure_path)
#         aligner = CEAligner()
#         aligner.set_reference(representative_structure)
#         target_name = f"{target['pdb_id']}_{target['chains']}_{target['model']}"
#         target_structure_path = os.path.join(pdb_chains_path, f"{target_name}.cif")
#         target_structure = parser.get_structure(target_name, target_structure_path)
#
#         aligner.align(target_structure)
#         rms = aligner.rms
#
#         result = {
#             'cluster_entry_id': target['id'],
#             'rms': rms
#         }
#
#         align_task_logger.info("Alignment completed successfully.")
#         return result
#     except Exception as e:
#         error_object = {'cluster_entry_id': target['id'], 'error_message': e}
#         return error_object
#

def map_representatives_to_targets(cluster_entries, pdb_chains_path):
    """
    Maps representative protein structures to their corresponding target structures for alignment.

    This function prepares a list of pairs consisting of representative and target protein structures. It is used to
    facilitate the structural alignment process, ensuring each target structure is aligned with the correct representative.

    Args:
        cluster_entries (list): A list of dictionaries containing cluster information for protein structures.
        pdb_chains_path (str): The file path where the protein chain files are stored.

    Returns:
        list: A list of dictionaries, each containing a pair of representative and target structures for alignment.
    """
    representative_map = {entry['cluster_id']: entry for entry in cluster_entries if entry['is_representative']}
    pairs_list = []

    for entry in cluster_entries:

        if not entry['is_representative']:
            print(representative_map)
            print(entry)
            rep = representative_map[entry['cluster_id']]
            pair_dict = {'representative': rep, 'target': entry, 'pdb_chains_path': pdb_chains_path}
            pairs_list.append(pair_dict)

    return pairs_list
