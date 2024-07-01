import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from http.client import HTTPException

import pandas as pd
import requests
from Bio import ExPASy, SwissProt
from rq import Queue
from sqlalchemy import func, exists
from urllib.parse import quote

from protein_metamorphisms_is.base.pipeline import PipelineBase
from protein_metamorphisms_is.helpers.logger.logger import setup_logger
from protein_metamorphisms_is.helpers.parser.parser import extract_float, process_chain_string
from protein_metamorphisms_is.information_system.base.extractor import ExtractorBase
from protein_metamorphisms_is.sql.model import Accession, ProteinGOTermAssociation, GOTerm, PDBReference, Sequence, \
    Protein


class UniProtExtractor(PipelineBase):
    """
    Handles the extraction and processing of data from UniProt. This includes retrieving protein sequences,
    annotations, and functional data. The class extends ExtractorBase to utilize base functionalities and ensures
    integration with a database to store and manage the extracted data.
    """

    def __init__(self, conf):
        """
        Initializes the extractor with configuration for logging and database access.
        Args:
            conf (dict): Configuration dictionary specifying operational parameters.
        """
        super().__init__(conf, session_required=True)
        self.logger.info("UniProtExtractor initialized with configuration.")

    def set_tasks(self):
        tasks = ['load_accesions', 'download_record', 'store_entry']
        workers = [1, 12, 1]

        for task,num_workers in zip(tasks, workers):
            self.queues[task] = Queue(task,connection=self.redis_conn)


    def set_targets(self):
        """
        Determines the data extraction targets from either a CSV file or a live API based on provided configuration.
        """
        csv_path = self.conf.get("load_accesion_csv")
        accession_column = self.conf.get("load_accesion_column")
        search_criteria = self.conf.get("search_criteria")
        limit = self.conf.get("limit")
        csv_tag = self.conf.get("csv_tag")

        if csv_path:
            self._load_access_from_csv(csv_path, accession_column, csv_tag)
        else:
            self.logger.warning("CSV path not provided; attempting to fetch accessions using API.")
            if search_criteria:
                self._fetch_accessions_from_api(search_criteria, limit)
            else:
                self.logger.error("No valid search criteria provided. Data extraction cannot proceed.")

    def fetch(self):
        """
        Initiates the download of UniProt entries using Redis queues to manage the download tasks.
        """
        self.logger.info("Starting the download of UniProt entries.")

        accessions = self.session.query(Accession).all()
        self.logger.info(f"Total proteins to download: {len(accessions)}")

        max_workers = self.conf.get("max_workers", 10)  # Default to 10 if not specified in config

        jobs = []
        for accession in accessions:
            job = self.process_queue.enqueue(download_record, accession.accession_code)
            jobs.append(job)

        while jobs:
            finished_jobs = [job for job in jobs if job.get_status() in ['finished', 'failed']]
            for job in finished_jobs:
                try:
                    data = job.result
                    if data:
                        self.data_queue.put(data)  # This should be changed if data_queue is not a multiprocessing.Queue
                        self.logger.info(f"Record for accession {job.args[0]} added to the queue.")
                    else:
                        self.logger.warning(f"No data found for accession {job.args[0]}")
                except Exception as e:
                    self.logger.error(f"Error processing the entry for accession {job.args[0]}: {e}")
                jobs.remove(job)

        # Signal to the consuming task(s) that fetching is complete by adding None or a similar sentinel value
        self.data_queue.put(None)
        self.logger.info("All data fetching and queuing completed.")
