import logging
import multiprocessing
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
from multiprocessing import pool, Pool

from Bio.PDB import CEAligner, MMCIFParser
from sqlalchemy.orm import sessionmaker
from tqdm import tqdm

from protein_data_handler.operations.base.bioinfo_operator import BioinfoOperatorBase
from protein_data_handler.sql.model import PDBChains, Cluster, PDBReference, CEAlignResults, CEAlignQueue


class CEAlign(BioinfoOperatorBase):
    """
    This class performs structural alignment of protein data using the Combinatorial Extension (CE) algorithm.
    CE algorithm is a popular method for protein structure alignment, known for its effectiveness in identifying
    resemblances between proteins that share very little sequence similarity. It works by creating an optimal
    alignment between two protein structures based on their backbone atom positions.
    """

    def __init__(self, conf):
        """
        Initializes the CEAlign object with configuration settings.

        Args:
            conf (dict): Configuration parameters, including database connections and operational settings.
        """
        super().__init__(conf, session_required=True)
        self.logger.info("CEAlign instance created with configuration.")

    def start(self):
        """
        Begins the structural alignment process. This method orchestrates the entire alignment workflow,
        starting from loading the cluster representatives, performing CE alignment, and handling any exceptions
        that occur during the process.
        """
        try:
            self.logger.info("Starting structural alignment process.")

            self.get_update_queue()

            while not self.check_empty_queue():
                batch_map = self.fetch_queue_items()
                alignment_map = map_representatives_to_targets(batch_map, self.conf.get('pdb_chains_path'))
                self.ce_align(alignment_map)


        except Exception as e:
            self.logger.error(f"Error during structural alignment process: {e}")
            raise

    def get_update_queue(self):
        """
        Update the alignment queue by adding new cluster entries and resetting stale entries.

        This method performs two main functions:
        1. Adds new clusters to the alignment queue (`CEAlignQueue`) that are not already present in it.
           It checks for clusters that are not in the queue and adds them with an initial 'pending' state.

        2. Resets the state of stale entries in the queue. Entries that are either in process or have
           encountered an error and haven't been updated within a specified timeout (default 60 minutes)
           are considered stale. Their state is reset to 'pending' to allow reprocessing.

        After processing, all changes are committed to the database.

        Note:
            - The 'state' values used in this method are 0 for 'pending', 1 for 'in process', 2 for 'completed',
              and 3 for 'error'.
        """

        self.logger.info("Loading cluster entries from database.")

        # Obtener una lista de IDs de cluster que ya están en la cola
        queued_cluster_ids = self.session.query(CEAlignQueue.cluster_entry_id).all()
        self.logger.info(f"Found {len(queued_cluster_ids)} clusters already in queue.")
        queued_cluster_ids = [c[0] for c in queued_cluster_ids]  # Convertir a una lista simple

        # Obtener clusters que no están en la cola
        clusters_not_queued = self.session.query(Cluster).filter(
            Cluster.id.notin_(queued_cluster_ids),
            Cluster.is_representative == False  # Add this condition
        ).all()

        self.logger.info(f"Found {len(clusters_not_queued)} clusters not in queue, adding to queue.")

        for cluster in clusters_not_queued:
            new_queue_item = CEAlignQueue(
                cluster_entry_id=cluster.id,
                state=0,  # Suponiendo que 0 es el estado para 'Pendiente'
                retry_count=0
            )
            self.session.add(new_queue_item)

        self.logger.info("New clusters added to queue. Checking for stale entries.")

        retry_timeout_date = datetime.now() - timedelta(minutes=self.conf.get('retry_timeout', 60))

        stale_entries = self.session.query(CEAlignQueue).filter(
            CEAlignQueue.state.in_([1, 3]),
            CEAlignQueue.updated_at < retry_timeout_date
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
        Check if the alignment queue is empty.

        Returns:
            bool: True if the queue is empty, False otherwise.
        """
        # Contar el número de elementos en la cola
        queue_count = self.session.query(CEAlignQueue).filter(
            CEAlignQueue.state.in_([0, 3]),
            CEAlignQueue.retry_count < self.conf.get('retry_count', 5)
        ).count()

        # Verificar si la cola está vacía
        is_empty = queue_count == 0

        if is_empty:
            self.logger.info("The alignment queue is empty.")
        else:
            self.logger.info(f"The alignment queue has {queue_count} items.")

        return is_empty

    def fetch_queue_items(self):
        batch_size = self.conf.get('batch_size', 100)  # Default to 100 if not specified

        # Fetch items from the queue
        queue_items = self.session.query(CEAlignQueue).filter(
            CEAlignQueue.state.in_([0, 3]),
            CEAlignQueue.retry_count < self.conf.get('retry_count', 5)
        ).order_by(
            CEAlignQueue.updated_at.desc()
        ).limit(
            batch_size
        ).all()

        self.logger.info(f"Fetched {len(queue_items)} items from the queue for processing.")

        # Update state to 1 (in progress) for each fetched item
        for item in queue_items:
            item.state = 1

        # Commit the state changes to the database
        self.session.commit()
        self.logger.info("Updated the state of fetched items to 'in progress'.")

        # Fetch related cluster information for the fetched queue items
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

    def ce_align(self, alignment_map):
        num_workers = self.conf.get('max_workers', 4)
        timeout = self.conf.get('task_timeout', 1)  # Un tiempo de espera de 300 segundos (5 minutos) por defecto
        self.logger.info(f"Performing CE alignment with {num_workers} workers.")

        results = []
        with Pool(num_workers) as pool:
            # Crear una lista de tareas asincrónicas
            tasks = [pool.apply_async(align_task, args=(entry,)) for entry in alignment_map]

            # Recoger los resultados con manejo de excepciones para timeout
            for task in tasks:
                try:
                    result = task.get(timeout=timeout)
                    results.append(result)
                except multiprocessing.TimeoutError:
                    self.logger.error("Timeout error: La tarea excedió el tiempo límite establecido")
                    # Agregar un manejo de error o resultado predeterminado si es necesario
                    # Por ejemplo: results.append({'error_message': 'Timeout'})

        self.logger.info(f"CE alignment completed for {len(results)} entries.")

        self.insert_ce_align_results(results)
        self.logger.info("CE alignment results inserted into the database successfully.")

    def insert_ce_align_results(self, results):
        """
        Inserts the results of the CE alignment into the database and updates the queue.

        Args:
            results (list): A list of dictionaries containing the alignment results or exceptions.
        """
        self.logger.info("Processing CE alignment results.")

        for result in results:
            if 'error_message' in result.keys():
                # Log error and update queue item to state 3 (error) and increment retry count
                self.logger.warn(f"Error in alignment task: {result}")

                # Update the queue item
                queue_item = self.session.query(CEAlignQueue).filter_by(
                    cluster_entry_id=result['cluster_entry_id']).first()
                if queue_item:
                    queue_item.state = 3  # Error state
                    queue_item.retry_count += 1
                    queue_item.error_message = str(result['error_message'])
                    self.session.add(queue_item)
            else:
                # Process successful result
                new_result = CEAlignResults(
                    cluster_entry_id=result['cluster_entry_id'],
                    rms=result['rms']
                )
                self.session.add(new_result)

                # Update the queue item to state 2 (completed)
                queue_item = self.session.query(CEAlignQueue).filter_by(
                    cluster_entry_id=result['cluster_entry_id']).first()
                if queue_item:
                    queue_item.state = 2  # Completed state
                    self.session.add(queue_item)

        # Commit the changes to the database
        self.session.commit()
        self.logger.info("CE alignment results processed and queue updated.")


def align_task(alignment_entry):
    align_task_logger = logging.getLogger("align_task")
    try:
        align_task_logger.info("Aligning structures...")
        representative = alignment_entry['representative']
        target = alignment_entry['target']
        pdb_chains_path = alignment_entry['pdb_chains_path']
        parser = MMCIFParser()
        representative_name = f"{representative['pdb_id']}_{representative['chains']}_{representative['model']}"
        representative_structure_path = os.path.join(pdb_chains_path, f"{representative_name}.cif")
        representative_structure = parser.get_structure(representative_name, representative_structure_path)
        aligner = CEAligner()
        aligner.set_reference(representative_structure)
        target_name = f"{target['pdb_id']}_{target['chains']}_{target['model']}"
        target_structure_path = os.path.join(pdb_chains_path, f"{target_name}.cif")
        target_structure = parser.get_structure(target_name, target_structure_path)

        print(target_name,representative_name)
        aligner.align(target_structure)
        rms = aligner.rms

        result = {
            'cluster_entry_id': target['id'],
            'rms': rms
        }

        align_task_logger.info("Alignment completed successfully.")
        return result
    except Exception as e:
        error_object = {'cluster_entry_id': target['id'], 'error_message': e}
        return error_object


def map_representatives_to_targets(cluster_entries, pdb_chains_path):
    representative_map = {entry['cluster_id']: entry for entry in cluster_entries if entry['is_representative']}
    pairs_list = []

    for entry in cluster_entries:

        if not entry['is_representative']:
            rep = representative_map[entry['cluster_id']]
            pair_dict = {'representative': rep, 'target': entry, 'pdb_chains_path': pdb_chains_path}
            pairs_list.append(pair_dict)

    return pairs_list
