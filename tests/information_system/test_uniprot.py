import unittest
from concurrent.futures import Future
from http.client import HTTPException
from unittest.mock import patch, MagicMock, call

import pandas as pd
import requests
from sqlalchemy import exists

from protein_metamorphisms_is.information_system.uniprot import UniProtExtractor
from protein_metamorphisms_is.sql.model import Accession, Protein, Sequence, PDBReference, GOTerm, \
    ProteinGOTermAssociation


class TestUniProtExtractor(unittest.TestCase):

    def setUp(self):
        self.conf = {
            "DB_USERNAME": "usuario",
            "DB_PASSWORD": "clave",
            "DB_HOST": "localhost",
            "DB_PORT": 5432,
            "DB_NAME": "BioData",
            "search_criteria": "(structure_3d:true)",
            "limit": 100,
            "max_workers": 2,
            "load_accesion_column": "accession",  # Asegúrate de que estos valores estén correctamente definidos
            "csv_tag": "test_tag"  # Asegúrate de que estos valores estén correctamente definidos
        }
        with patch('protein_metamorphisms_is.information_system.base.extractor.setup_logger') as mock_logger_setup, \
                patch('protein_metamorphisms_is.sql.base.database_manager.create_engine'), \
                patch('protein_metamorphisms_is.sql.base.database_manager.sessionmaker',
                      return_value=MagicMock()):
            mock_logger_setup.return_value = MagicMock()  # Configuramos un MagicMock como logger
            self.extractor = UniProtExtractor(self.conf)
            self.extractor.logger = MagicMock()  # Asegura que el logger se inicializa adecuadamente
            self.extractor.session = MagicMock()
            self.extractor.data_queue = MagicMock()

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._load_access_from_csv')
    def test_set_targets_with_csv(self, mock_load_csv):
        # Configura la prueba para usar una ruta de CSV
        self.extractor.conf['load_accesion_csv'] = "path/to/csv"
        self.extractor.set_targets()
        mock_load_csv.assert_called_once_with("path/to/csv", "accession", "test_tag")

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._fetch_accessions_from_api')
    def test_set_targets_with_api(self, mock_fetch_api):
        # No se proporciona la ruta del CSV, debería intentar cargar desde la API
        self.extractor.set_targets()
        mock_fetch_api.assert_called_once_with("(structure_3d:true)", 100)

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._fetch_accessions_from_api')
    def test_set_targets_no_valid_search_criteria(self, mock_fetch_api):
        self.extractor.conf['search_criteria'] = None
        self.extractor.set_targets()

    @patch('pandas.read_csv')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._process_new_accessions')
    def test_load_access_from_csv_success(self, mock_process_new_accessions, mock_read_csv):
        # Simula que read_csv retorna un DataFrame específico
        mock_data = pd.DataFrame({
            'accession': ['A12345', 'B67890', 'C13579', 'A12345']  # Datos duplicados para probar 'unique'
        })
        mock_read_csv.return_value = mock_data

        csv_path = "path/to/csv"
        accession_column = "accession"
        csv_tag = "test_tag"

        # Llama a la función con los parámetros mockeados
        self.extractor._load_access_from_csv(csv_path, accession_column, csv_tag)
        # # Verifica que se lea el CSV correctamente
        mock_read_csv.assert_called_once_with(csv_path)

        self.extractor.logger.info.assert_called_with("Loaded 3 unique accession codes from CSV.")

    @patch('pandas.read_csv')
    def test_load_access_from_csv_failure(self, mock_read_csv):
        # Simula una excepción al intentar leer el CSV
        mock_read_csv.side_effect = Exception("Failed to read CSV")

        csv_path = "path/to/invalid_csv"
        accession_column = "accession"
        csv_tag = "test_tag"

        # Llama a la función y espera una excepción
        self.extractor._load_access_from_csv(csv_path, accession_column, csv_tag)

        # Verifica que el log de error sea correcto
        self.extractor.logger.error.assert_called_once()

        # Asegúrate de que la excepción capturada sea logueada adecuadamente
        error_message = self.extractor.logger.error.call_args[0][0]
        self.assertIn("Failed to load or process CSV", error_message)

    @patch('requests.get')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._process_new_accessions')
    def test_fetch_accessions_from_api_success(self, mock_process_new_accessions, mock_requests_get):
        # Simula una respuesta exitosa de la API
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "A12345\nB67890\nC13579"
        mock_requests_get.return_value = mock_response

        search_criteria = "(structure_3d:true)"
        limit = 100

        # Llama a la función con los parámetros mockeados
        self.extractor._fetch_accessions_from_api(search_criteria, limit)

        # Verifica que se haga la solicitud HTTP correctamente
        encoded_search_criteria = requests.utils.quote(search_criteria)
        expected_url = f"https://rest.uniprot.org/uniprotkb/stream?query={encoded_search_criteria}&format=list&size={limit}"
        mock_requests_get.assert_called_once_with(expected_url)

        # Verifica que se procesen los códigos de acceso obtenidos de la API
        mock_process_new_accessions.assert_called_once_with(['A12345', 'B67890', 'C13579'],
                                                            self.extractor.conf.get("tag"))

        # Verifica que el log de información sea correcto
        self.extractor.logger.info.assert_any_call(f"Fetching data from URL: {expected_url}")
        self.extractor.logger.info.assert_any_call("Retrieved 3 accessions from UniProt API.")

    @patch('requests.get')
    def test_fetch_accessions_from_api_failure(self, mock_requests_get):
        # Simula una excepción al hacer la solicitud HTTP
        mock_requests_get.side_effect = requests.RequestException("Failed to fetch data")

        search_criteria = "(structure_3d:true)"
        limit = 100

        # Llama a la función y espera una excepción
        self.extractor._fetch_accessions_from_api(search_criteria, limit)

        # Verifica que el log de error sea correcto
        self.extractor.logger.error.assert_called_once()

        # Asegúrate de que la excepción capturada sea logueada adecuadamente
        error_message = self.extractor.logger.error.call_args[0][0]
        self.assertIn("Failed to fetch data from UniProt", error_message)

    def test_process_new_accessions(self):
        accessions = ['A12345', 'B67890', 'C13579']
        tag = "test_tag"

        # Simula los códigos de acceso existentes en la base de datos
        existing_accessions_query = MagicMock()
        existing_accessions_query.filter.return_value.all.return_value = [('A12345',)]
        self.extractor.session.query.return_value = existing_accessions_query

        # Llama a la función con los parámetros mockeados
        self.extractor._process_new_accessions(accessions, tag)

        # Verifica que se haya registrado el log de información correcto
        self.extractor.logger.info.assert_any_call("Processing 3 accessions.")

        # Verifica que se haya consultado la base de datos correctamente
        self.extractor.session.query.assert_called_once_with(Accession.accession_code)
        existing_accessions_query.filter.assert_called_once()

        # Verifica que se hayan guardado los nuevos códigos de acceso
        new_accessions = [Accession(accession_code='B67890', primary=True, tag=tag),
                          Accession(accession_code='C13579', primary=True, tag=tag)]
        self.extractor.session.bulk_save_objects.assert_called_once()
        # # Verifica que se haya hecho commit en la sesión
        self.extractor.session.commit.assert_called_once()
        #
        # # Verifica que se haya registrado el log de información de nuevos accesos
        self.extractor.logger.info.assert_any_call("Added 2 new accessions to the database.")

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._download_record')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor.store_entry')
    def test_fetch(self, mock_store_entry, mock_download_record):
        mock_download_record.side_effect = [MagicMock(), None, Exception("Failed to download")]

        # Configura la respuesta de la base de datos simulada con tres accesiones
        self.extractor.session.query.return_value.all.return_value = [
            MagicMock(accession_code="ABC123"),
            MagicMock(accession_code="DEF456"),
            MagicMock(accession_code="GHI789"),
        ]

        self.extractor.fetch()
        # Verifica que el logger registró los mensajes de inicio y finalización correctamente
        self.extractor.logger.info.assert_any_call("Starting the download of UniProt entries.")
        self.extractor.logger.info.assert_any_call("Total proteins to download: 3")
        self.extractor.logger.info.assert_any_call("Record for accession ABC123 added to the queue.")
        self.extractor.logger.warning.assert_any_call("No data found for accession DEF456")
        self.extractor.logger.error.assert_any_call(
            "Error processing the entry for accession GHI789: Failed to download")
        self.extractor.logger.info.assert_any_call("All data fetching and queuing completed.")

        # Verifica que el método _download_record se llamó con los códigos de acceso correctos
        expected_calls = [call("ABC123"), call("DEF456"), call("GHI789")]
        mock_download_record.assert_has_calls(expected_calls, any_order=True)

    @patch('protein_metamorphisms_is.information_system.uniprot.ExPASy.get_sprot_raw')
    @patch('protein_metamorphisms_is.information_system.uniprot.SwissProt.read')
    def test_download_record_success(self, mock_swissprot_read, mock_get_sprot_raw):
        mock_handle = MagicMock()
        mock_record = MagicMock()
        mock_get_sprot_raw.return_value = mock_handle
        mock_swissprot_read.return_value = mock_record

        result = self.extractor._download_record('ABC123')

        mock_get_sprot_raw.assert_called_once_with('ABC123')
        mock_swissprot_read.assert_called_once_with(mock_handle)
        self.assertEqual(result, mock_record)

    @patch('protein_metamorphisms_is.information_system.uniprot.ExPASy.get_sprot_raw')
    @patch('protein_metamorphisms_is.information_system.uniprot.SwissProt.read')
    def test_download_record_failure(self, mock_swissprot_read, mock_get_sprot_raw):
        mock_get_sprot_raw.side_effect = Exception('Mock Exception')

        result = self.extractor._download_record('ABC123')

        mock_get_sprot_raw.assert_called_once_with('ABC123')
        mock_swissprot_read.assert_not_called()
        self.extractor.logger.error.assert_called_once_with('Error downloading the entry ABC123: Mock Exception')
        self.assertIsNone(result)

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._get_or_create_protein')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._update_protein_details')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._handle_accessions')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._handle_cross_references')
    def test_store_entry_success(self, mock_handle_cross_references, mock_handle_accessions,
                                 mock_update_protein_details, mock_get_or_create_protein):
        mock_data = MagicMock()
        mock_data.entry_name = "ABC123"
        mock_protein = MagicMock()
        mock_get_or_create_protein.return_value = mock_protein

        self.extractor.store_entry(mock_data)

        self.extractor.logger.debug.assert_called_once_with(
            "Attempting to store or update protein data for entry_name ABC123.")
        mock_get_or_create_protein.assert_called_once_with(mock_data)
        mock_update_protein_details.assert_called_once_with(mock_protein, mock_data)
        mock_handle_accessions.assert_called_once_with(mock_protein, mock_data.accessions)
        mock_handle_cross_references.assert_called_once_with(mock_protein, mock_data.cross_references)
        self.extractor.session.commit.assert_called_once()
        self.extractor.logger.info.assert_called_once_with("Successfully stored or updated data for protein ABC123.")

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._get_or_create_protein')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._update_protein_details')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._handle_accessions')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._handle_cross_references')
    def test_store_entry_http_exception(self, mock_handle_cross_references, mock_handle_accessions,
                                        mock_update_protein_details, mock_get_or_create_protein):
        mock_data = MagicMock()
        mock_data.entry_name = "ABC123"
        self.extractor._get_or_create_protein.side_effect = HTTPException("Mock HTTP Exception")

        self.extractor.store_entry(mock_data)

        self.extractor.logger.error.assert_called_once_with(
            "HTTP error occurred while processing ABC123: Mock HTTP Exception")
        self.extractor.session.rollback.assert_called_once()

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._get_or_create_protein')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._update_protein_details')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._handle_accessions')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._handle_cross_references')
    def test_store_entry_general_exception(self, mock_handle_cross_references, mock_handle_accessions,
                                           mock_update_protein_details, mock_get_or_create_protein):
        mock_data = MagicMock()
        mock_data.entry_name = "ABC123"
        self.extractor._get_or_create_protein.side_effect = Exception("Mock General Exception")

        self.extractor.store_entry(mock_data)

        self.extractor.logger.error.assert_called_once_with(
            "An unexpected error occurred while processing ABC123: Mock General Exception")
        self.extractor.session.rollback.assert_called_once()

    def test_get_or_create_protein_existing(self):
        mock_data = MagicMock()
        mock_data.entry_name = "ABC123"
        mock_protein = MagicMock()
        self.extractor.session.query.return_value.filter_by.return_value.first.return_value = mock_protein

        result = self.extractor._get_or_create_protein(mock_data)

        self.extractor.session.query.return_value.filter_by.assert_called_once_with(entry_name="ABC123")
        self.extractor.session.query.return_value.filter_by.return_value.first.assert_called_once()
        self.assertEqual(result, mock_protein)
        self.extractor.session.add.assert_not_called()

    def test_get_or_create_protein_new(self):
        mock_data = MagicMock()
        mock_data.entry_name = "ABC123"
        self.extractor.session.query.return_value.filter_by.return_value.first.return_value = None

        result = self.extractor._get_or_create_protein(mock_data)

        self.extractor.session.query.return_value.filter_by.assert_called_once_with(entry_name="ABC123")
        self.extractor.session.query.return_value.filter_by.return_value.first.assert_called_once()
        self.assertIsInstance(result, Protein)
        self.assertEqual(result.entry_name, "ABC123")
        self.extractor.session.add.assert_called_once_with(result)

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._get_or_create_sequence')
    def test_update_protein_details(self, mock_get_or_create_sequence):
        mock_sequence = MagicMock()
        mock_get_or_create_sequence.return_value = mock_sequence

        mock_protein = MagicMock()
        mock_data = MagicMock()
        mock_data.sequence = 'mock_sequence'
        mock_data.data_class = 'mock_data_class'
        mock_data.molecule_type = 'mock_molecule_type'
        mock_data.sequence_length = 123
        mock_data.created = [MagicMock()]
        mock_data.sequence_update = [MagicMock()]
        mock_data.annotation_update = [MagicMock()]
        mock_data.description = 'mock_description'
        mock_data.gene_name = 'mock_gene_name'
        mock_data.organism = 'mock_organism'
        mock_data.organelle = 'mock_organelle'
        mock_data.organism_classification = ['class1', 'class2']
        mock_data.taxonomy_id = ['id1', 'id2']
        mock_data.host_organism = ['host1', 'host2']
        mock_data.host_taxonomy_id = ['host_id1', 'host_id2']
        mock_data.comments = ['comment1', 'comment2']
        mock_data.keywords = 'mock_keywords'
        mock_data.protein_existence = 1
        mock_data.seqinfo = 'mock_seqinfo'

        self.extractor._update_protein_details(mock_protein, mock_data)

        mock_get_or_create_sequence.assert_called_once_with('mock_sequence')
        self.assertEqual(mock_protein.sequence, mock_sequence)
        self.assertEqual(mock_protein.data_class, 'mock_data_class')
        self.assertEqual(mock_protein.molecule_type, 'mock_molecule_type')
        self.assertEqual(mock_protein.sequence_length, 123)
        self.assertEqual(mock_protein.created_date, mock_data.created[0])
        self.assertEqual(mock_protein.sequence_update_date, mock_data.sequence_update[0])
        self.assertEqual(mock_protein.annotation_update_date, mock_data.annotation_update[0])
        self.assertEqual(mock_protein.description, 'mock_description')
        self.assertEqual(mock_protein.gene_name, 'mock_gene_name')
        self.assertEqual(mock_protein.organism, 'mock_organism')
        self.assertEqual(mock_protein.organelle, 'mock_organelle')
        self.assertEqual(mock_protein.organism_classification, 'class1,class2')
        self.assertEqual(mock_protein.taxonomy_id, 'id1,id2')
        self.assertEqual(mock_protein.host_organism, 'host1,host2')
        self.assertEqual(mock_protein.host_taxonomy_id, 'host_id1,host_id2')
        self.assertEqual(mock_protein.comments, 'comment1; comment2')
        self.assertEqual(mock_protein.keywords, 'mock_keywords')
        self.assertEqual(mock_protein.protein_existence, 1)
        self.assertEqual(mock_protein.seqinfo, 'mock_seqinfo')
        self.extractor.session.add.assert_called_once_with(mock_protein)

    def test_get_or_create_sequence_existing(self):
        mock_sequence = "MSEQ123"
        mock_existing_sequence = MagicMock()
        self.extractor.session.query.return_value.filter_by.return_value.first.return_value = mock_existing_sequence

        result = self.extractor._get_or_create_sequence(mock_sequence)

        self.extractor.session.query.return_value.filter_by.assert_called_once_with(sequence=mock_sequence)
        self.extractor.session.query.return_value.filter_by.return_value.first.assert_called_once()
        self.assertEqual(result, mock_existing_sequence)
        self.extractor.session.add.assert_not_called()

    def test_get_or_create_sequence_new(self):
        mock_sequence = "MSEQ123"
        self.extractor.session.query.return_value.filter_by.return_value.first.return_value = None

        result = self.extractor._get_or_create_sequence(mock_sequence)

        self.extractor.session.query.return_value.filter_by.assert_called_once_with(sequence=mock_sequence)
        self.extractor.session.query.return_value.filter_by.return_value.first.assert_called_once()
        self.assertIsInstance(result, Sequence)
        self.assertEqual(result.sequence, mock_sequence)
        self.extractor.session.add.assert_called_once_with(result)

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._get_or_create_accession')
    def test_handle_accessions(self, mock_get_or_create_accession):
        mock_protein = MagicMock()
        mock_protein.entry_name = "PROT123"
        mock_accessions = ["ACC123", "ACC456"]

        mock_accession1 = MagicMock()
        mock_accession2 = MagicMock()
        mock_get_or_create_accession.side_effect = [mock_accession1, mock_accession2]

        self.extractor._handle_accessions(mock_protein, mock_accessions)

        mock_get_or_create_accession.assert_has_calls([call("ACC123"), call("ACC456")])
        self.assertEqual(mock_accession1.protein_entry_name, "PROT123")
        self.assertEqual(mock_accession2.protein_entry_name, "PROT123")
        self.extractor.session.add.assert_has_calls([call(mock_accession1), call(mock_accession2)])

    def test_get_or_create_accession_existing(self):
        mock_accession_code = "ACC123"
        mock_existing_accession = MagicMock()

        self.extractor.session.query.return_value.scalar.return_value = True
        self.extractor.session.query.return_value.filter.return_value.first.return_value = mock_existing_accession

        result = self.extractor._get_or_create_accession(mock_accession_code)

        self.assertEqual(result, mock_existing_accession)
        self.extractor.session.add.assert_not_called()  # Verifica que no se llama a add para un acceso existente

    def test_get_or_create_accession_new(self):
        mock_accession_code = "ACC123"

        self.extractor.session.query.return_value.scalar.return_value = False

        result = self.extractor._get_or_create_accession(mock_accession_code)

        self.assertIsInstance(result, Accession)
        self.assertFalse(result.primary)
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._handle_pdb_reference')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._handle_go_reference')
    def test_handle_cross_references(self, mock_handle_go_reference, mock_handle_pdb_reference):
        mock_protein = MagicMock()
        mock_cross_references = [
            ("PDB", "PDB123"),
            ("GO", "GO456"),
            ("OTHER", "OTHER789")
        ]

        self.extractor._handle_cross_references(mock_protein, mock_cross_references)

        mock_handle_pdb_reference.assert_called_once_with(mock_protein, ("PDB", "PDB123"))
        mock_handle_go_reference.assert_called_once_with(mock_protein, ("GO", "GO456"))
        mock_handle_pdb_reference.assert_called_with(mock_protein, ("PDB", "PDB123"))
        mock_handle_go_reference.assert_called_with(mock_protein, ("GO", "GO456"))
        self.assertEqual(mock_handle_pdb_reference.call_count, 1)
        self.assertEqual(mock_handle_go_reference.call_count, 1)

    @patch('protein_metamorphisms_is.information_system.uniprot.process_chain_string')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._get_or_create_pdb_reference')
    def test_handle_pdb_reference(self, mock_get_or_create_pdb_reference, mock_process_chain_string):
        mock_protein = MagicMock()
        mock_protein.sequence.sequence = "SEQUENCE"
        mock_reference = ["PDB", "some_id", "some_method", "1.5", "A:1-7"]

        mock_process_chain_string.return_value = ("A", 1, 7)
        mock_pdb_ref = MagicMock()
        mock_get_or_create_pdb_reference.return_value = mock_pdb_ref

        self.extractor._handle_pdb_reference(mock_protein, mock_reference)

        mock_process_chain_string.assert_called_once_with("A:1-7")
        mock_get_or_create_pdb_reference.assert_called_once_with(mock_reference, "SEQUENC")
        self.assertEqual(mock_pdb_ref.protein, mock_protein)
        self.extractor.session.add.assert_called_once_with(mock_pdb_ref)

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._get_or_create_sequence')
    @patch('protein_metamorphisms_is.information_system.uniprot.extract_float')
    def test_get_or_create_pdb_reference_new(self, mock_extract_float, mock_get_or_create_sequence):
        mock_reference = ["PDB", "PDB123", "X-ray", "1.5"]
        mock_sequence = "SEQUENCE"
        mock_new_sequence = MagicMock()
        mock_get_or_create_sequence.return_value = mock_new_sequence
        mock_extract_float.return_value = 1.5

        self.extractor.session.query.return_value.scalar.return_value = False

        result = self.extractor._get_or_create_pdb_reference(mock_reference, mock_sequence)

        # self.extractor.session.query.assert_called_once_with(exists().where(PDBReference.pdb_id == "PDB123"))
        mock_get_or_create_sequence.assert_called_once_with(mock_sequence)
        mock_extract_float.assert_called_once_with("1.5")
        self.assertIsInstance(result, PDBReference)
        self.assertEqual(result.pdb_id, "PDB123")
        self.assertEqual(result.method, "X-ray")
        self.assertEqual(result.resolution, 1.5)
        self.assertEqual(result.sequence, mock_new_sequence)

    def test_get_or_create_pdb_reference_existing(self):
        mock_reference = ["PDB", "PDB123", "X-ray", "1.5"]
        mock_sequence = "SEQUENCE"
        mock_existing_pdb_ref = MagicMock()

        self.extractor.session.query.return_value.scalar.return_value = True
        self.extractor.session.query.return_value.filter.return_value.first.return_value = mock_existing_pdb_ref

        result = self.extractor._get_or_create_pdb_reference(mock_reference, mock_sequence)

        # self.extractor.session.query.assert_called_once_with(exists().where(PDBReference.pdb_id == "PDB123"))
        # self.extractor.session.query.return_value.filter.assert_called_once_with(PDBReference.pdb_id == "PDB123")
        self.assertEqual(result, mock_existing_pdb_ref)
        self.extractor.session.add.assert_not_called()

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._get_or_create_go_term')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor._get_or_create_association')
    def test_handle_go_reference(self, mock_get_or_create_association, mock_get_or_create_go_term):
        mock_protein = MagicMock()
        mock_protein.entry_name = "PROT123"
        mock_reference = ["GO", "GO:0008150", "biological_process", "IEA"]

        mock_go_term = MagicMock()
        mock_go_term.go_id = "GO:0008150"
        mock_get_or_create_go_term.return_value = mock_go_term

        mock_association = None
        mock_get_or_create_association.return_value = mock_association

        with patch.object(self.extractor.logger, 'info') as mock_logger_info:
            self.extractor._handle_go_reference(mock_protein, mock_reference)

            mock_get_or_create_go_term.assert_called_once_with(mock_reference)
            mock_get_or_create_association.assert_called_once_with(mock_protein.entry_name, mock_go_term.go_id)
            mock_logger_info.assert_called_once_with(
                f"Association between {mock_protein.entry_name} and GO Term {mock_go_term.go_id} already exists.")

    def test_get_or_create_go_term_existing(self):
        mock_reference = ["GO", "GO:0008150", "biological_process:description", "IEA:other"]
        mock_go_term = MagicMock()

        self.extractor.session.query.return_value.filter_by.return_value.first.return_value = mock_go_term

        result = self.extractor._get_or_create_go_term(mock_reference)

        self.extractor.session.query.return_value.filter_by.assert_called_once_with(go_id="GO:0008150")
        self.assertEqual(result, mock_go_term)
        self.extractor.session.add.assert_not_called()

    def test_get_or_create_go_term_new(self):
        mock_reference = ["GO", "GO:0008150", "biological_process:description", "IEA:other"]

        self.extractor.session.query.return_value.filter_by.return_value.first.return_value = None

        result = self.extractor._get_or_create_go_term(mock_reference)

        self.extractor.session.query.return_value.filter_by.assert_called_once_with(go_id="GO:0008150")
        self.assertIsInstance(result, GOTerm)
        self.assertEqual(result.go_id, "GO:0008150")
        self.assertEqual(result.category, "biological_process")
        self.assertEqual(result.description, "description")
        self.assertEqual(result.evidences, "IEA")
        self.extractor.session.add.assert_called_once_with(result)

    def test_get_or_create_association_existing(self):
        entry_name = "PROT123"
        go_id = "GO:0008150"
        mock_existing_association = MagicMock()

        self.extractor.session.query.return_value.scalar.return_value = True
        self.extractor.session.query.return_value.filter.return_value.first.return_value = mock_existing_association

        result = self.extractor._get_or_create_association(entry_name, go_id)

        self.extractor.session.add.assert_not_called()

    def test_get_or_create_association_new(self):
        entry_name = "PROT123"
        go_id = "GO:0008150"

        self.extractor.session.query.return_value.scalar.return_value = False

        result = self.extractor._get_or_create_association(entry_name, go_id)

        self.assertIsInstance(result, ProteinGOTermAssociation)
        self.assertEqual(result.protein_entry_name, entry_name)
        self.assertEqual(result.go_id, go_id)
        self.extractor.session.add.assert_called_once_with(result)