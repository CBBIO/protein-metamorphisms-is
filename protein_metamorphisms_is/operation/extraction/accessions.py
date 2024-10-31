import traceback
from urllib.parse import quote
import pandas as pd
import requests

from protein_metamorphisms_is.sql.model.entities.protein.accesion import Accession
from protein_metamorphisms_is.tasks.base import BaseTaskInitializer


class AccessionManager(BaseTaskInitializer):
    """
    The AccessionManager class is responsible for managing accession data within the system.
    It extends the BaseTaskInitializer to handle tasks such as loading accession data from
    CSV files and fetching accession codes from the UniProt API.

    **Purpose**

    The AccessionManager class provides functionalities to manage and process accession codes
    for biological data, ensuring they are properly stored and maintained in the database.

    **Key Features**

    - **Load from CSV**: Load accession data from a specified CSV file and process it.
    - **Fetch from API**: Fetch accession data from the UniProt API based on specified search criteria.
    - **Database Integration**: Seamlessly integrates with the database to store and manage accession codes.
    - **Logging**: Inherits logging capabilities from BaseTaskInitializer for tracking operations and errors.

    **Customization**

    Subclasses can override the enqueue, process, and store_entry methods to implement
    specific task logic related to accession data management.

    Attributes:
        conf (dict): Configuration dictionary loaded from YAML or other sources.
        logger (Logger): Logger instance for logging task-specific information.
        session (Session): Database session used for ORM operations.

    Example Usage:

    .. code-block:: python

       from protein_metamorphisms_is.tasks.accessions import AccessionManager

       class MyAccessionManager(AccessionManager):
           def enqueue(self):
               # Implementation of enqueue logic
               pass

           def process(self, _):
               # Implementation of process logic
               pass

           def store_entry(self, record):
               # Implementation of store_entry logic
               pass
    """

    def __init__(self, conf, session_required=True):
        """
        Initialize the AccessionManager.

        This constructor sets up the logger, configuration, and database session (if required).

        Args:
            conf (dict): Configuration dictionary.
            session_required (bool): Whether a database session is required.
                                     If True, the session is initialized.
        """
        super().__init__(conf, session_required)

    def load_accessions_from_csv(self):
        """
        Loads accessions from a specified CSV file and processes them for data fetching.

        This method reads accession codes from a CSV file, ensures they are unique,
        and processes them by invoking the `_process_new_accessions` method.

        Raises:
            Exception: If there is an issue loading or processing the CSV file.
        """
        try:
            csv_path = self.conf['load_accesion_csv']
            accession_column = self.conf['load_accesion_column']
            csv_tag = self.conf['tag']

            data = pd.read_csv(csv_path)
            accessions = data[accession_column].dropna().unique()[:1000
                         ]
            self.logger.info(f"Loaded {len(accessions)} unique accession codes from CSV.")
            self._process_new_accessions(accessions, csv_tag)
        except Exception as e:
            self.logger.error(f"Failed to load or process CSV: {traceback.format_exc()}")

    def fetch_accessions_from_api(self):
        """
        Fetches accession codes from the UniProt API based on the specified search criteria.

        This method sends a request to the UniProt API using the configured search criteria
        and processes the resulting accession codes by invoking the `_process_new_accessions` method.

        Raises:
            requests.RequestException: If there is an issue with the API request.
        """
        try:
            search_criteria = self.conf['search_criteria']
            limit = self.conf['limit']
            tag = self.conf.get('tag')

            encoded_search_criteria = quote(search_criteria)
            url = f"https://rest.uniprot.org/uniprotkb/stream?query={encoded_search_criteria}&format=list&size={limit}"
            self.logger.info(f"Fetching data from URL: {url}")

            response = requests.get(url)
            response.raise_for_status()
            accessions = response.text.strip().split("\n")[:3000]
            self.logger.info(f"Retrieved {len(accessions)} accessions from UniProt API.")
            self._process_new_accessions(accessions, tag)
        except requests.RequestException as e:
            self.logger.error(f"Failed to fetch data from UniProt: {e}")

    def _process_new_accessions(self, accessions, tag):
        """
        Processes newly fetched accession codes, checking against the database to avoid duplicates,
        and saves them if they are new.

        Args:
            accessions (list): List of accession codes to process.
            tag (str): Tag to associate with new accessions.

        Raises:
            Exception: If there is an issue with database operations.
        """
        self.logger.info(f"Processing {len(accessions)} accessions.")
        existing_accessions = {acc[0] for acc in self.session.query(Accession.accession_code).filter(
            Accession.accession_code.in_(accessions)).all()}
        new_accessions = [Accession(accession_code=acc, primary=True, tag=tag) for acc in accessions if
                          acc not in existing_accessions]
        self.session.bulk_save_objects(new_accessions)
        self.session.commit()
        self.logger.info(f"Added {len(new_accessions)} new accessions to the database.")

    def enqueue(self):
        """
        Abstract method for enqueuing tasks.

        This method should be implemented in subclasses to define how accession tasks are enqueued.
        """
        pass

    def store_entry(self, record):
        """
        Abstract method for storing processed entries.

        This method should be implemented in subclasses to define how processed accession entries
        are stored in the database.
        """
        pass

    def process(self, _):
        """
        Abstract method for processing tasks.

        This method should be implemented in subclasses to define the logic for processing
        accession data.
        """
        pass

    def start(self):
        """
        Abstract method for starting the task processing.

        This method should be implemented in subclasses to define how the task processing
        is initiated for accession management.
        """
        pass
