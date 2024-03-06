from unittest import mock, TestCase
from unittest.mock import patch, MagicMock, call
import unittest

from sqlalchemy.orm import Session

# Asegúrate de que la importación de PDBExtractor sea correcta según tu estructura de paquetes
from protein_metamorphisms_is.information_system.pdb import PDBExtractor
from protein_metamorphisms_is.information_system.pdb import ChainSelect


class TestPDBExtractor(unittest.TestCase):
    def setUp(self):
        self.conf = {
            "DB_USERNAME": "usuario",
            "DB_PASSWORD": "clave",
            "DB_HOST": "localhost",
            "DB_PORT": 5432,
            "DB_NAME": "BioData",
            "search_criteria": "(structure_3d:true)",
            "limit": 100,
            "max_workers": 2  # Asegúrate de definir esto si tu clase lo requiere
        }
        with patch('protein_metamorphisms_is.information_system.base.extractor.setup_logger') as mock_logger, \
                patch('protein_metamorphisms_is.sql.base.database_manager.create_engine'), \
                patch('protein_metamorphisms_is.sql.base.database_manager.sessionmaker',
                      return_value=MagicMock()):
            self.mock_logger = mock_logger
            self.pdb_id = "sample_pdb_id"
            self.pdb_reference = MagicMock(pdb_id=self.pdb_id)
            self.extractor = PDBExtractor(self.conf)
            self.pdb_reference_id = 12345
            self.local_session = mock.MagicMock(spec=Session)
            self.pdb_file_path = "path/to/test/file.cif"

    @patch('protein_metamorphisms_is.sql.model.PDBReference')
    @patch('protein_metamorphisms_is.information_system.pdb.or_')
    def test_load_pdb_ids(self, mock_or, mock_pdb_reference):
        # Setup
        expected_resolution_threshold = 2.0
        expected_pdb_references = [MagicMock(), MagicMock()]  # Mock PDBReference objects

        self.extractor.conf = {"resolution_threshold": expected_resolution_threshold}
        self.extractor.session = MagicMock()
        self.extractor.session.query.return_value.filter.return_value.all.return_value = expected_pdb_references

        # Mock the or_ function to simply pass its arguments through
        mock_or.side_effect = lambda *args: args

        # Execute
        result = self.extractor.load_pdb_ids()

        # Verify
        mock_or.assert_called()
        self.assertEqual(result, expected_pdb_references)

    @patch('protein_metamorphisms_is.information_system.pdb.PDBExtractor.download_pdb_structures')
    @patch('protein_metamorphisms_is.information_system.pdb.PDBExtractor.load_pdb_ids')
    def test_start_success(self, mock_load_pdb_ids, mock_download_pdb_structures):
        # Setup
        mock_pdb_references = [MagicMock(), MagicMock()]
        mock_load_pdb_ids.return_value = mock_pdb_references

        # Execute
        self.extractor.start()

        # Verify
        mock_load_pdb_ids.assert_called_once()
        mock_download_pdb_structures.assert_called_once()

    @patch('protein_metamorphisms_is.information_system.pdb.PDBExtractor.download_pdb_structures')
    @patch('protein_metamorphisms_is.information_system.pdb.PDBExtractor.load_pdb_ids')
    def test_start_exception(self, mock_load_pdb_ids, mock_download_pdb_structures):
        # Setup
        mock_load_pdb_ids.side_effect = Exception("Test exception")

        # Execute and Verify

        self.extractor.start()

        # The logger's error method should be called once with the exception message
        self.extractor.logger.error.assert_called()

    @patch('protein_metamorphisms_is.information_system.pdb.sessionmaker')  # Adjust the import path as necessary
    @patch('Bio.PDB.PDBList')
    @patch('protein_metamorphisms_is.information_system.pdb.PDBExtractor.populate_pdb_chains')
    def test_download_and_process_pdb_structure(self, mock_populate_pdb_chains, mock_pdb_list, mock_sessionmaker):
        # Setup
        mock_session = MagicMock()
        mock_sessionmaker.return_value = mock_session

        mock_pdb_list_instance = MagicMock()
        mock_pdb_list.return_value = mock_pdb_list_instance
        mock_pdb_list_instance.retrieve_pdb_file.return_value = self.pdb_file_path

        # Execute
        self.extractor.download_and_process_pdb_structure(self.pdb_reference)

        # Verify
        mock_sessionmaker.assert_called_once_with(bind=self.extractor.engine)
        mock_pdb_list.assert_called_once_with(server="ftp.wwpdb.org", pdb="pdb_files")
        mock_pdb_list_instance.retrieve_pdb_file.assert_called_once_with(self.pdb_id, file_format="mmCif",
                                                                         pdir="pdb_files")
        mock_populate_pdb_chains.assert_called()

    @patch('protein_metamorphisms_is.information_system.pdb.sessionmaker')  # Adjust the import path as necessary
    @patch('Bio.PDB.PDBList')
    @patch('protein_metamorphisms_is.information_system.pdb.PDBExtractor.populate_pdb_chains')
    def test_download_and_process_pdb_structure_with_exception(self, mock_populate_pdb_chains, mock_pdb_list,
                                                               mock_sessionmaker):
        # Setup
        mock_session = MagicMock()
        mock_sessionmaker.return_value = mock_session

        mock_pdb_list_instance = MagicMock()
        mock_pdb_list.return_value = mock_pdb_list_instance
        mock_pdb_list_instance.retrieve_pdb_file.side_effect = Exception("Download failed")

        # Execute and Verify

        self.extractor.download_and_process_pdb_structure(self.pdb_reference)

        # Verify that an exception leads to the correct logging
        # This assumes you have a logger configured in your PDBExtractor class
        with patch.object(self.extractor.logger, 'error') as mock_log_error:
            self.extractor.download_and_process_pdb_structure(self.pdb_reference)
            mock_log_error.assert_called()

    @patch('protein_metamorphisms_is.information_system.pdb.MMCIFParser')
    @patch('os.path.isfile', return_value=False)
    @patch('protein_metamorphisms_is.information_system.pdb.MMCIFIO')
    @patch('protein_metamorphisms_is.information_system.pdb.PDBChains', autospec=True)
    def test_populate_pdb_chains(self, mock_pdb_chains, mock_mmcifio, mock_isfile, mock_mmcif_parser):
        # Setup mock structure with correct iterable behavior
        mock_structure = MagicMock()
        mock_chain = MagicMock()
        mock_chain.get_id.return_value = 'A'
        mock_model = MagicMock()
        mock_model.get_id.return_value = 0  # Ensure model ID is set to 0

        # Correctly simulate residues in the chain for sequence generation
        mock_residue_ala = MagicMock()
        mock_residue_ala.id = (' ', 1, ' ')
        mock_residue_ala.resname = 'ALA'
        mock_residue_gly = MagicMock()
        mock_residue_gly.id = (' ', 2, ' ')
        mock_residue_gly.resname = 'GLY'

        # Make the structure and model iterable
        mock_structure.__iter__.return_value = iter([mock_model])
        mock_model.__iter__.return_value = iter([mock_chain])
        mock_chain.__iter__.return_value = iter([mock_residue_ala,
                                                 mock_residue_gly])

        mock_parser_instance = MagicMock()
        mock_parser_instance.get_structure.return_value = mock_structure
        mock_mmcif_parser.return_value = mock_parser_instance

        # Assume pdb_reference_id exists in the database
        self.local_session.query.return_value.filter.return_value.first.return_value = (self.pdb_reference_id,)

        # Execute
        self.extractor.populate_pdb_chains(self.pdb_file_path, self.pdb_reference_id, self.local_session)

        # Verify
        mock_mmcif_parser.assert_called_once_with()
        mock_parser_instance.get_structure.assert_called_once_with(self.pdb_reference_id, self.pdb_file_path)
        self.local_session.add.assert_called()
        mock_mmcifio.assert_called_once()
        mock_isfile.assert_called()

        # Ensure the mock call matches the expected call, including the sequence
        mock_pdb_chains.assert_called_with(chains='A', sequence='AG', pdb_reference_id=self.pdb_reference_id, model=0)

        # Verify file saving
        mock_mmcifio_instance = mock_mmcifio.return_value
        mock_mmcifio_instance.set_structure.assert_called_once_with(mock_structure)
        mock_mmcifio_instance.save.assert_called()

        # Verify database operations
        self.local_session.commit.assert_called_once()
        self.local_session.close.assert_called_once()

    @patch('protein_metamorphisms_is.information_system.pdb.MMCIFParser')
    @patch('os.path.isfile', return_value=False)
    @patch('protein_metamorphisms_is.information_system.pdb.MMCIFIO')
    @patch('protein_metamorphisms_is.information_system.pdb.PDBChains', autospec=True)
    def test_populate_pdb_chains_with_no_pdb_reference(self, mock_pdb_chains, mock_mmcifio, mock_isfile,
                                                       mock_mmcif_parser):
        # Setup mock structure with correct iterable behavior
        mock_structure = MagicMock()
        mock_chain = MagicMock()
        mock_chain.get_id.return_value = 'A'
        mock_model = MagicMock()
        mock_model.get_id.return_value = 0  # Ensure model ID is set to 0

        # Correctly simulate residues in the chain for sequence generation
        mock_residue_ala = MagicMock()
        mock_residue_ala.id = (' ', 1, ' ')
        mock_residue_ala.resname = 'ALA'
        mock_residue_gly = MagicMock()
        mock_residue_gly.id = (' ', 2, ' ')
        mock_residue_gly.resname = 'GLY'

        # Make the structure and model iterable
        mock_structure.__iter__.return_value = iter([mock_model])
        mock_model.__iter__.return_value = iter([mock_chain])
        mock_chain.__iter__.return_value = iter([mock_residue_ala, mock_residue_gly])

        mock_parser_instance = MagicMock()
        mock_parser_instance.get_structure.return_value = mock_structure
        mock_mmcif_parser.return_value = mock_parser_instance

        # Simulate the case where the pdb_reference_id does not exist in the database
        self.local_session.query.return_value.filter.return_value.first.return_value = None

        # Execute
        self.extractor.populate_pdb_chains(self.pdb_file_path, self.pdb_reference_id, self.local_session)

        # Verify
        mock_mmcif_parser.assert_called_once_with()
        mock_parser_instance.get_structure.assert_called_once_with(self.pdb_reference_id, self.pdb_file_path)
        self.local_session.add.assert_called()
        mock_mmcifio.assert_called_once()
        mock_isfile.assert_called()

        # Ensure the mock call matches the expected call, including the sequence
        # This time, pdb_reference_id_value should be None in the call
        mock_pdb_chains.assert_called_with(chains='A', sequence='AG', pdb_reference_id=None, model=0)

        # Verify file saving
        mock_mmcifio_instance = mock_mmcifio.return_value
        mock_mmcifio_instance.set_structure.assert_called_once_with(mock_structure)
        mock_mmcifio_instance.save.assert_called()

        # Verify database operations
        self.local_session.commit.assert_called_once()
        self.local_session.close.assert_called_once()

    @mock.patch('protein_metamorphisms_is.information_system.pdb.ThreadPoolExecutor')
    @mock.patch('protein_metamorphisms_is.information_system.pdb.PDBExtractor.download_and_process_pdb_structure')
    def test_download_pdb_structures(self, mock_download_and_process, mock_executor):
        mock_executor.return_value.__enter__.return_value.map = mock.MagicMock()
        pdb_references = [mock.MagicMock(), mock.MagicMock()]
        self.extractor.download_pdb_structures(pdb_references)

        # Verify ThreadPoolExecutor is called with the correct number of workers
        mock_executor.assert_called_once_with(max_workers=self.conf["max_workers"])

        # Verify download_and_process_pdb_structure is scheduled for execution for each PDB reference
        args_list = mock_executor.return_value.__enter__.return_value.map.call_args[0]
        self.assertEqual(args_list[0], self.extractor.download_and_process_pdb_structure)
        self.assertEqual(list(args_list[1]), pdb_references)

        # Verify logging
        self.extractor.logger.info.assert_called_with(
            f"Downloading PDB structures with {self.conf['max_workers']} workers")


class TestChainSelect(TestCase):
    def setUp(self):
        self.chain_id = 'A'
        self.model_id = 0
        self.chain_select = ChainSelect(chain_id=self.chain_id, model_id=self.model_id)

    def test_accept_chain_true(self):
        mock_chain = mock.MagicMock()
        mock_chain.get_id.return_value = self.chain_id

        self.assertTrue(self.chain_select.accept_chain(mock_chain))

    def test_accept_chain_false(self):
        mock_chain = mock.MagicMock()
        mock_chain.get_id.return_value = 'B'


        self.assertFalse(self.chain_select.accept_chain(mock_chain))

    def test_accept_model_true(self):
        mock_model = mock.MagicMock()
        mock_model.get_id.return_value = self.model_id

        self.assertTrue(self.chain_select.accept_model(mock_model))

    def test_accept_model_false(self):
        mock_model = mock.MagicMock()
        mock_model.get_id.return_value = 1

        self.assertFalse(self.chain_select.accept_model(mock_model))
