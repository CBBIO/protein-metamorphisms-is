import gzip
import traceback
from io import BytesIO
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

    **Configuration**

    :param conf: Configuration dictionary loaded from YAML or other sources.
    :type conf: dict
    :param session_required: Whether a database session is required. Default is True.
    :type session_required: bool

    Example Usage:
        Below are examples showing how to use the key functionalities of the `AccessionManager`:

        .. code-block:: python

            from protein_metamorphisms_is.tasks.accessions import AccessionManager

            # Initialize the AccessionManager with configuration
            config = {
                'load_accesion_csv': 'path_to_csv_file.csv',
                'load_accesion_column': 'accession_column_name',
                'tag': 'example_tag',
                'search_criteria': 'example_criteria',
                'limit': 200,
                'debug': True
            }
            accession_manager = AccessionManager(config)

            # Load accessions from CSV
            accession_manager.load_accessions_from_csv()

            # Fetch accessions from API
            accession_manager.fetch_accessions_from_api()

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
            limit_execution = self.conf.get('limit_execution', None)  # Obtener el límite de ejecución
            data = pd.read_csv(csv_path)
            accessions = list(data[accession_column].dropna().unique())

            # Aplica el límite si está configurado
            if limit_execution and isinstance(limit_execution, int):
                accessions = accessions[:limit_execution]

            self.logger.info(f"Loaded {len(accessions)} unique accession codes from CSV.")
            self._process_new_accessions(accessions, csv_tag)
        except Exception:
            self.logger.error(f"Failed to load or process CSV: {traceback.format_exc()}")

    def fetch_accessions_from_api(self):
        """
        Fetches accession codes from the UniProt API based on the specified search criteria.

        This method queries the UniProt API with compression and pagination enabled. It processes
        the results and stores new accession codes in the database.

        :raises requests.RequestException: If there is an error in the API request.
        :raises gzip.BadGzipFile: If the response is not a valid gzip file.
        :raises Exception: For any other errors during processing.
        """
        try:
            search_criteria = self.conf['search_criteria']
            limit_per_page = self.conf.get('limit', 100)  # Límite por página
            total_limit = self.conf.get('limit_execution', None)  # Límite total de ejecución
            tag = self.conf.get('tag')

            encoded_search_criteria = quote(search_criteria)
            base_url = f"https://rest.uniprot.org/uniprotkb/search?compressed=true&format=list&query={encoded_search_criteria}&size={limit_per_page}"

            all_accessions = []
            next_cursor = None

            while True:
                # Construye la URL con el cursor si está disponible
                url = base_url if not next_cursor else f"{base_url}&cursor={next_cursor}"
                self.logger.info(f"Fetching data from URL: {url}")

                response = requests.get(url)
                response.raise_for_status()

                # Descomprime el contenido usando gzip
                with gzip.GzipFile(fileobj=BytesIO(response.content)) as f:
                    decompressed_content = f.read().decode('utf-8')

                # Obtén los accesiones y añádelos a la lista total
                accessions = decompressed_content.strip().split("\n")

                all_accessions.extend(accessions)

                # Verifica si alcanzamos el límite total
                if total_limit and len(all_accessions) >= total_limit:
                    all_accessions = all_accessions[:total_limit]
                    break

                # Revisa si hay más páginas para obtener el siguiente cursor
                link_header = response.headers.get("link")
                if link_header and 'rel="next"' in link_header:
                    next_cursor = link_header.split("cursor=")[-1].split(">")[0]
                else:
                    break  # No hay más páginas, termina el bucle
            self.logger.info(f"Total accessions retrieved: {len(all_accessions)}")
            self._process_new_accessions(all_accessions, tag)
        except requests.RequestException as e:
            self.logger.error(f"Failed to fetch data from UniProt: {e}")
        except gzip.BadGzipFile:
            self.logger.error("Received an invalid gzip file. Unable to decompress the data.")
        except Exception as e:
            self.logger.error(f"An error occurred while processing the response: {e}")

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
        existing_accessions = {acc[0] for acc in self.session.query(Accession.code).filter(
            Accession.code.in_(accessions)).all()}
        new_accessions = [Accession(code=acc, primary=True, tag=tag) for acc in accessions if
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
        (Not used)
        """
        pass

    def process(self, _):
        """
        (Not used)
        """
        pass

    def start(self):
        """
        (Not used)
        """
        pass
