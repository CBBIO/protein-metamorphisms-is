import logging
import multiprocessing
import re
import subprocess

import os
from datetime import datetime, timedelta
from multiprocessing import Pool

from Bio.PDB import CEAligner, MMCIFParser

from protein_data_handler.operations.base.operator import OperatorBase
from protein_data_handler.sql.model import PDBChains, Cluster, PDBReference, TMAlignResults, TMAlignQueue


class TMAlign(OperatorBase):
    """
    """

    def __init__(self, conf):
        """
        Initializes the TMAlign object with configuration settings.

        Args:
            conf (dict): Configuration parameters, including database connections and operational settings.
        """
        super().__init__(conf, session_required=True)
        self.logger.info("TMAlign instance created with configuration.")

    def start(self):
        """
        Begin the structural alignment process using the CE algorithm.

        This method manages the workflow of the alignment process, including loading clusters, executing alignments,
        and handling any exceptions encountered during the process. Progress and errors are logged appropriately.
        """
        try:
            self.logger.info("Starting TM-align structural alignment process.")

            self.get_update_queue()

            while not self.check_empty_queue():
                batch_map = self.fetch_queue_items()
                alignment_map = map_representatives_to_targets(batch_map, self.conf.get('pdb_chains_path'))
                self.tm_align(alignment_map)

        except Exception as e:
            self.logger.error(f"Error during structural alignment process: {e}")
            raise

    def get_update_queue(self):
        """
        Updates the alignment queue by adding new cluster entries and resetting stale entries.

        This method performs two primary functions in managing the CE alignment queue (`CEAlignQueue`):
        1. Addition of New Clusters: Identifies clusters not currently in the queue and adds them with an initial
           'pending' state (state 0). This ensures that all relevant clusters are queued for alignment processing.

        2. Resetting Stale Entries: Updates the status of stale entries in the queue. An entry is considered stale if
           it is either 'in process' (state 1) or 'error' (state 3) and has not been updated within a specified timeout
           (default is 60 minutes). These stale entries are reset to 'pending' (state 0) for reprocessing.

        The method concludes by committing all changes to the database, thereby keeping the queue up-to-date.

        Note:
            - State values in the `CEAlignQueue` are defined as 0 for 'pending', 1 for 'in process', 2 for 'completed',
              and 3 for 'error'.
        """
        self.logger.info("Loading cluster entries from database.")

        queued_cluster_ids = self.session.query(TMAlignQueue.cluster_entry_id).all()
        self.logger.info(f"Found {len(queued_cluster_ids)} clusters already in queue.")
        queued_cluster_ids = [c[0] for c in queued_cluster_ids]

        clusters_not_queued = self.session.query(Cluster).filter(
            Cluster.id.notin_(queued_cluster_ids),
            Cluster.is_representative == False  # Add this condition
        ).all()

        self.logger.info(f"Found {len(clusters_not_queued)} clusters not in queue, adding to queue.")

        for cluster in clusters_not_queued:
            new_queue_item = TMAlignQueue(
                cluster_entry_id=cluster.id,
                state=0,
                retry_count=0
            )
            self.session.add(new_queue_item)

        self.logger.info("New clusters added to queue. Checking for stale entries.")

        retry_timeout_date = datetime.now() - timedelta(minutes=self.conf.get('retry_timeout', 60))

        stale_entries = self.session.query(TMAlignQueue).filter(
            TMAlignQueue.state.in_([1, 3]),
            TMAlignQueue.updated_at < retry_timeout_date
        ).all()
        self.logger.info(f"Found {len(stale_entries)} stale queue entries, resetting state to 'pending'.")

        for entry in stale_entries:
            entry.state = 0
            entry.retry_count += 1
            self.session.add(entry)

        self.session.commit()
        self.logger.info("Queue update process completed.")

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
        queue_count = self.session.query(TMAlignQueue).filter(
            TMAlignQueue.state.in_([0, 3]),
            TMAlignQueue.retry_count < self.conf.get('retry_count', 5)
        ).count()

        is_empty = queue_count == 0

        if is_empty:
            self.logger.info("The alignment queue is empty.")
        else:
            self.logger.info(f"The alignment queue has {queue_count} items.")

        return is_empty

    def fetch_queue_items(self):
        """
        Fetches a batch of items from the alignment queue for processing.

        This method retrieves a set number of items from the CE alignment queue (`CEAlignQueue`). The items fetched are
        either in a 'pending' (state 0) or 'error' (state 3) state and have not exceeded the maximum retry count. The
        method updates the state of these items to 'in progress' (state 1) and commits these changes to the database.

        Additionally, it retrieves relevant cluster information for each queue item. This includes fetching both the
        representative and non-representative cluster entries, which are essential for the subsequent alignment process.

        Returns:
            list: A list of dictionaries, each containing detailed information about a cluster entry, including
                  its ID, cluster ID, representativeness, chain identifiers, model numbers, and PDB IDs.
        """
        batch_size = self.conf.get('batch_size', 100)

        queue_items = self.session.query(TMAlignQueue).filter(
            TMAlignQueue.state.in_([0, 3]),
            TMAlignQueue.retry_count < self.conf.get('retry_count', 5)
        ).order_by(
            TMAlignQueue.updated_at.desc()
        ).limit(
            batch_size
        ).all()

        self.logger.info(f"Fetched {len(queue_items)} items from the queue for processing.")

        for item in queue_items:
            item.state = 1

        self.session.commit()
        self.logger.info("Updated the state of fetched items to 'in progress'.")

        cluster_ids = [item.cluster_entry_id for item in queue_items]
        clusters = self.session.query(
            Cluster.id,
            Cluster.cluster_id,
            Cluster.is_representative,
            PDBChains.chains,
            PDBChains.model,
            PDBReference.pdb_id
        ).join(
            PDBChains, Cluster.pdb_chain_id == PDBChains.id
        ).join(
            PDBReference, PDBChains.pdb_reference_id == PDBReference.id
        ).filter(
            Cluster.id.in_(cluster_ids)
        ).all()

        self.logger.info(f"Fetched cluster information for {len(clusters)} queue items.")

        # Fetch representative clusters corresponding to the non-representative ones
        representative_cluster_ids = [cluster.cluster_id for cluster in clusters if not cluster.is_representative]
        representative_clusters = self.session.query(
            Cluster.id,
            Cluster.cluster_id,
            Cluster.is_representative,
            PDBChains.chains,
            PDBChains.model,
            PDBReference.pdb_id
        ).join(
            PDBChains, Cluster.pdb_chain_id == PDBChains.id
        ).join(
            PDBReference, PDBChains.pdb_reference_id == PDBReference.id
        ).filter(
            Cluster.cluster_id.in_(representative_cluster_ids),
            Cluster.is_representative == True
        ).all()

        self.logger.info(f"Fetched {len(representative_clusters)} representative cluster entry.")

        # Combine both lists
        clusters.extend(representative_clusters)

        self.logger.info(f"Total loaded cluster information for {len(clusters)} items.")

        # Prepare results
        results = []
        for row in clusters:
            row_dict = {
                "id": row[0],
                "cluster_id": row[1],
                "is_representative": row[2],
                "chains": row[3],
                "model": row[4],
                "pdb_id": row[5]
            }
            results.append(row_dict)

        return results

    def tm_align(self, alignment_map):
        """
        Executes structural alignments using the TM-align algorithm, in a parallelized manner.

        This method leverages multiprocessing to perform structural alignments concurrently, improving efficiency
        for large datasets. It distributes alignment tasks across multiple worker processes and aggregates the
        results. If a task exceeds the specified timeout, an error is logged, and the task is handled appropriately.

        Args:
            alignment_map (list): A list of dictionaries, each containing the details required for an alignment task.

        The method:
        1. Initializes a pool of worker processes for parallel execution of alignment tasks.
        2. Submits alignment tasks asynchronously and waits for their completion, adhering to a specified timeout.
        3. Collects results from each completed task and handles any exceptions, such as timeouts.
        4. Inserts the alignment results into the database for further analysis or reference.

        Note:
            - The number of worker processes ('max_workers') and the task timeout ('task_timeout') are configurable.
            - In case of a timeout, an appropriate error message is logged, and the task is processed accordingly.
        """
        num_workers = self.conf.get('max_workers', 4)
        timeout = self.conf.get('task_timeout', 1)
        self.logger.info(f"Performing TM-align with {num_workers} workers.")

        results = []
        with Pool(num_workers) as pool:
            tasks = [pool.apply_async(align_task, args=(entry,)) for entry in alignment_map]

            for task in tasks:
                try:
                    result = task.get(timeout=timeout)
                    results.append(result)
                except multiprocessing.TimeoutError:
                    self.logger.error("Timeout error: La tarea excedió el tiempo límite establecido")

        self.logger.info(f"TM-align completed for {len(results)} entries.")

        self.insert_tm_align_results(results)
        self.logger.info("TM-align results inserted into the database successfully.")

    def insert_tm_align_results(self, results):
        """
        Inserts the results of TM-align structural alignment tasks into the database and updates the alignment queue.

        This method processes the results obtained from the TM-align alignment tasks. For each result, it either stores the
        alignment data in the database or handles any errors that occurred during the alignment process. The method also
        updates the state of each queue item based on the outcome of its corresponding alignment task.

        Args:
            results (list): A list of dictionaries, each containing either the alignment results or error information.
                            The dictionaries include the cluster entry ID and, depending on success or failure,
                            the RMSD value, sequence identity, and TM-scores for both chains, or an error message.

        Process:
            - For successful alignments, the method stores the alignment results in the `TMAlignResults` table.
            - For failed alignments, indicated by an 'error_message' key in the result, the method logs the error and
              updates the queue item's state to 'error' (state 3), incrementing the retry count.
            - The method then updates the queue item's state to 'completed' (state 2) for successful alignments.

        Note:
            - The state of each queue item in `TMAlignQueue` is updated to reflect the result of the alignment task.
            - The states are defined as 0 for 'pending', 1 for 'in process', 2 for 'completed', and 3 for 'error'.
        """
        self.logger.info("Processing TM-align alignment results.")

        for result in results:
            if result['error_message']:
                self.logger.warn(f"Error in alignment task: {result['error_message']}")

                queue_item = self.session.query(TMAlignQueue).filter_by(
                    cluster_entry_id=result['cluster_entry_id']).first()
                if queue_item:
                    queue_item.state = 3
                    queue_item.retry_count += 1
                    queue_item.error_message = result['error_message']
                    self.session.add(queue_item)
            else:
                new_result = TMAlignResults(
                    cluster_entry_id=result['cluster_entry_id'],
                    rmsd=result['rmsd'],
                    seq_id=result['seq_id'],
                    tm_score_chain_1=result['tm_score_chain_1'],
                    tm_score_chain_2=result['tm_score_chain_2']
                )
                self.session.add(new_result)

                queue_item = self.session.query(TMAlignQueue).filter_by(
                    cluster_entry_id=result['cluster_entry_id']).first()
                if queue_item:
                    queue_item.state = 2  # Mark as completed
                    self.session.add(queue_item)

        self.session.commit()
        self.logger.info("TM-align alignment results processed and queue updated.")


def align_task(alignment_entry):
    """
    Performs the alignment task for a given pair of protein structures using TM-align.

    This function aligns a target protein structure with a representative structure using TM-align. It logs the
    process and handles any exceptions that occur during the alignment. The function is designed to be run asynchronously
    in a multiprocessing environment.

    Args:
        alignment_entry (dict): A dictionary containing the representative and target protein structures' information and paths.

    Returns:
        dict: A dictionary with the alignment result or error message, if any.
    """
    align_task_logger = logging.getLogger("align_task")
    try:
        align_task_logger.info("Aligning structures using TM-align...")
        representative = alignment_entry['representative']
        target = alignment_entry['target']
        pdb_chains_path = alignment_entry['pdb_chains_path']
        representative_name = f"{representative['pdb_id']}_{representative['chains']}_{representative['model']}"
        representative_structure_path = os.path.join(pdb_chains_path, f"{representative_name}.cif")
        target_name = f"{target['pdb_id']}_{target['chains']}_{target['model']}"
        target_structure_path = os.path.join(pdb_chains_path, f"{target_name}.cif")

        # Construye el comando para ejecutar TMalign
        command = ["TMalign", representative_structure_path, target_structure_path]

        # Ejecuta el comando
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Lee la salida y el error (si existe)
        stdout, stderr = process.communicate()

        if process.returncode == 0:
            # Procesa la salida si el comando fue exitoso
            # Expresiones regulares para buscar los valores deseados
            rmsd_pattern = re.compile(r"RMSD=\s*(\d+\.\d+)")
            seq_id_pattern = re.compile(r"Seq_ID=n_identical/n_aligned=\s*(\d+\.\d+)")
            tm_score_chain_1_pattern = re.compile(r"TM-score=\s*(\d+\.\d+)\s*\(if normalized by length of Chain_1\)")
            tm_score_chain_2_pattern = re.compile(r"TM-score=\s*(\d+\.\d+)\s*\(if normalized by length of Chain_2\)")

            rmsd = float(rmsd_pattern.search(stdout).group(1)) if rmsd_pattern.search(stdout) else None
            seq_id = float(seq_id_pattern.search(stdout).group(1)) if seq_id_pattern.search(stdout) else None
            tm_score_chain_1 = float(
                tm_score_chain_1_pattern.search(stdout).group(1)) if tm_score_chain_1_pattern.search(stdout) else None
            tm_score_chain_2 = float(
                tm_score_chain_2_pattern.search(stdout).group(1)) if tm_score_chain_2_pattern.search(stdout) else None

            result = {
                'cluster_entry_id': target['id'],
                'rmsd': rmsd,
                'seq_id': seq_id,
                'tm_score_chain_1': tm_score_chain_1,
                'tm_score_chain_2': tm_score_chain_2,
                'error_message': None
            }

        else:
            result = {
                'cluster_entry_id': target['id'],
                'error_message': stderr
            }

        align_task_logger.info("Alignment completed successfully.")
        return result
    except Exception as e:
        align_task_logger.error(f"Error during alignment task: {str(e)}")
        return {
            'cluster_entry_id': target['id'],
            'error_message': str(e)
        }


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
            rep = representative_map[entry['cluster_id']]
            pair_dict = {'representative': rep, 'target': entry, 'pdb_chains_path': pdb_chains_path}
            pairs_list.append(pair_dict)

    return pairs_list
