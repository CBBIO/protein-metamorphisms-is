import importlib
import multiprocessing
from datetime import datetime, timedelta
from multiprocessing import Pool

from sqlalchemy.orm import aliased

from protein_metamorphisms_is.operations.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import PDBChains, Cluster, PDBReference, StructuralAlignmentQueue, \
    StructuralAlignmentType, StructuralAlignmentResults


class StructuralAlignmentManager(OperatorBase):
    """
    Manages the structural alignment process of proteins using various alignment algorithms.

    This class extends `OperatorBase` to implement a task management system for structural alignment, allowing
    the execution of different alignment tasks that are configurable through dynamic modules.

    Attributes:
        conf (dict): Configuration of the instance, including database connections and operational settings.
    """

    def __init__(self, conf):
        """
        Initializes an instance of `StructuralAlignmentManager` with configuration settings.

        Args:
            conf (dict): Configuration parameters, including database connections and operational settings.
        """
        super().__init__(conf)
        self.logger.info("Structural Alignment Manager instance created.")

    def fetch_tasks_info(self):
        """
        Fetches and prepares alignment task modules based on the configuration.

        This method dynamically imports alignment task modules specified in the configuration and stores
        references to these modules in a dictionary for later use in the alignment process.
        """
        structural_alignment_types = self.session.query(StructuralAlignmentType).all()
        self.types = {}
        base_module_path = 'protein_metamorphisms_is.operations.structural_alignment_tasks'

        for type_obj in structural_alignment_types:
            if type_obj.id in self.conf['structural_alignment']['types']:
                # Construye el nombre completo del módulo
                module_name = f"{base_module_path}.{type_obj.task_name}"
                # Importa dinámicamente el módulo usando importlib
                module = importlib.import_module(module_name)
                # Almacena la referencia al módulo en el diccionario self.types
                self.types[type_obj.id] = module

    def start(self):
        """
        Begin the structural alignment process.

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
        """
        Updates the alignment queue with new tasks and manages stale entries.

        This method loads clusters from the database, updates the queue with new clusters not yet queued for alignment,
        and resets the state of stale entries in the queue.
        """
        for type_entry in self.types.items():
            type_id = type_entry[0]
            type_task_name = type_entry[1]
            self.logger.info("Loading cluster entries from database.")

            queued_cluster_ids = self.session.query(StructuralAlignmentQueue.cluster_entry_id) \
                .filter(StructuralAlignmentQueue.alignment_type_id == type_id) \
                .all()
            self.logger.info(f"Found a total of {len(queued_cluster_ids)} clusters already in queue.")
            queued_cluster_ids = [c[0] for c in queued_cluster_ids]

            clusters_not_queued = self.session.query(Cluster).filter(
                Cluster.id.notin_(queued_cluster_ids),
                not Cluster.is_representative
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
         Determines if the alignment queue is empty, indicating no pending or error tasks awaiting processing.

         Evaluates the alignment queue to assess the presence of tasks in 'pending' or 'error' status that have not
         yet exceeded the retry limit. This is crucial for deciding whether to fetch new tasks or proceed with retrying
         error tasks. The method specifically looks for tasks that are either awaiting initial processing (state 0)
         or have encountered errors but are eligible for retry (state 3), based on a retry count that is below
         a predefined threshold set in the configuration.

         Returns:
             bool: True if the queue is empty, indicating there are no tasks in 'pending' or 'error' status
                   that require processing. False otherwise, suggesting that there are tasks that need attention.
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
        """
         Retrieves a batch of alignment tasks from the queue, prioritizing those in 'pending' state.

         This method performs a complex query to the database to fetch a specified number of alignment tasks, up to
         the 'batch_size' limit defined in the configuration. It includes tasks that are either pending execution
         (state 0) or have previously encountered errors and are eligible for retry (state 3), provided they have
         not exceeded the maximum retry count also specified in the configuration.

         The method constructs a subquery to identify the representative structure for each cluster involved in the
         alignment tasks, facilitating a direct comparison between target structures and their respective representatives.
         It then joins this subquery with the main queue items query to ensure each task fetched includes comprehensive
         information about the target structure and its representative counterpart.

         Returns:
             list: A list of SQLAlchemy model instances, each representing a queue item. These instances include
                   detailed information about the task, such as the alignment type, cluster ID, PDB IDs, chain identifiers,
                   and model numbers for both the target structure and its representative.
         """
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
            Cluster.is_representative
        ).subquery()

        # Alias para hacer más fácil la referencia a la subconsulta
        rep_structures = aliased(representative_structures_subquery, name='rep_structures')

        # Realiza la consulta principal con un join a la subconsulta
        queue_items = self.session.query(
            StructuralAlignmentQueue.id.label("queue_entry_id"),
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
        """
        Executes alignment tasks for a batch of queue items using multiprocessing.

        This method leverages Python's multiprocessing capabilities to perform structural alignments in parallel,
        according to the number of workers specified in the configuration (`max_workers`). Each alignment task is
        executed asynchronously with a timeout limit (`task_timeout`), also specified in the configuration. The method
        handles both successful completions and timeouts by logging the outcomes and storing the results or errors.

        The alignment tasks are dynamically determined based on the alignment type associated with each queue item,
        allowing for flexibility in the alignment process. Results from completed tasks are collected and later
        inserted into the database.

        Args:
            queue_items (list): A list of queue items to be aligned. Each item contains necessary information
                                for performing the alignment, including the type of alignment to execute.

        Notes:
            - The method logs the start and completion of the alignment process, including any errors encountered.
            - Upon completion or timeout, results are passed to `insert_results` for database insertion.
        """
        num_workers = self.conf.get('max_workers', 4)
        timeout = self.conf['structural_alignment']['task_timeout']
        self.logger.info(f"Performing CE alignment with {num_workers} workers.")

        results = []

        with Pool(num_workers) as pool:
            tasks_with_args = []
            for item in queue_items:
                args = (item, self.conf)
                task = pool.apply_async(self.types[item.alignment_type_id].align_task, args=args)
                tasks_with_args.append((task, args))

            for task, args in tasks_with_args:
                try:
                    result = task.get(timeout=timeout)
                    results.append(result)
                except multiprocessing.TimeoutError:
                    error_message = f"Timeout error: La tarea excedió el tiempo límite establecido ({timeout})"
                    self.logger.error(error_message)
                    results.append((args[0].queue_entry_id, {'error_message': error_message}))

        self.logger.info(f"CE alignment completed for {len(results)} entries.")
        self.insert_results(results)
        self.logger.info("CE alignment results inserted into the database successfully.")

    def insert_results(self, results):
        """
        Inserts the outcomes of structural alignment tasks into the database and updates the status of queue items.

        After structural alignment tasks are completed, this method is responsible for processing the results,
        which may include successful alignment data or error messages for tasks that failed. Depending on the outcome,
        it updates the database with the alignment results for successful tasks or logs and records errors for tasks
        that encountered issues. Additionally, it updates the status of each task in the alignment queue to reflect
        its current state, whether completed successfully, failed with an error, or pending retry based on the retry
        policy defined in the configuration.

        Args:
            results (list): A list containing the results of the alignment tasks. Each element in the list is a tuple
                            where the first element is the queue entry ID, and the second element is a dictionary with
                            either the alignment results (e.g., RMS values) for successful tasks or an error message for
                            tasks that failed.

        Process:
            - Iterates through the list of results, processing each based on its content (success or failure).
            - For successful tasks, stores alignment results (e.g., RMS values) in the database.
            - For tasks that failed, logs the error and updates the task's status in the queue to 'error' (state 3),
              while also incrementing the retry count for possible future attempts.
            - Updates the task's status in the queue to 'completed' (state 2) for successful alignments.
            - Commits all changes to the database once all results are processed.

        Note:
            The method ensures that the alignment queue is accurately updated to reflect the outcome of each task,
            facilitating efficient management of pending, in-process, and completed tasks.
        """
        self.logger.info("Processing CE alignment results.")
        for result in results:
            queue_entry_id = result[0]
            result = result[1]
            if 'error_message' in result.keys():
                self.logger.warning(f"Error in alignment task: {result}")

                queue_item = self.session.query(StructuralAlignmentQueue).filter_by(
                    id=queue_entry_id).first()
                if queue_item:
                    queue_item.state = 3  # Error state
                    queue_item.retry_count += 1
                    queue_item.error_message = str(result['error_message'])
                    self.session.add(queue_item)
            else:
                new_result = StructuralAlignmentResults(
                    **result
                )
                self.session.add(new_result)

                queue_item = self.session.query(StructuralAlignmentQueue).filter_by(
                    id=queue_entry_id).first()
                if queue_item:
                    queue_item.state = 2
                    self.session.add(queue_item)

        self.session.commit()
        self.logger.info("Structural alignment results processed and queue updated.")
