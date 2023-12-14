import unittest
from unittest.mock import patch, MagicMock
from protein_data_handler.fasta import FastaHandler

from requests.exceptions import RequestException
import pandas as pd


class TestFastaHandler(unittest.TestCase):

    @patch('os.makedirs')
    @patch('os.path.exists')
    def test_init(self, mock_exists, mock_makedirs):
        mock_exists.return_value = False
        session_mock = MagicMock()
        data_dir = './test/data'
        output_dir = './test/output'
        handler = FastaHandler(session_mock, data_dir, output_dir)
        mock_makedirs.assert_called()
        self.assertEqual(handler.session, session_mock)
        self.assertEqual(handler.data_dir, data_dir)
        self.assertEqual(handler.output_dir, output_dir)

    @patch('os.path.exists')
    @patch('concurrent.futures.ThreadPoolExecutor')
    @patch('protein_data_handler.fasta.FastaHandler.download_and_store_fasta')
    def test_download_fastas(self, mock_download_and_store_fasta, mock_executor, mock_exists):
        mock_exists.return_value = False
        mock_download_and_store_fasta.return_value = ('path', [('PDB1', 'chain1', 'id1', 'sequence1')])
        session_mock = MagicMock()
        handler = FastaHandler(session_mock, './data', './output')
        pdb_ids = ['PDB1', 'PDB2']
        handler.download_fastas(pdb_ids)
        self.assertEqual(mock_download_and_store_fasta.call_count, len(pdb_ids))

    @patch('os.path.exists')
    @patch('concurrent.futures.ThreadPoolExecutor')
    @patch('protein_data_handler.fasta.FastaHandler.download_and_store_fasta')
    def test_download_fastas_invalid_input(self, mock_download_and_store_fasta, mock_executor, mock_exists):
        mock_exists.return_value = False
        session_mock = MagicMock()
        handler = FastaHandler(session_mock, './data', './output')

        # Prueba con un argumento no válido (no es una lista)
        with self.assertRaises(ValueError):
            handler.download_fastas("not-a-list")

        # Prueba con una lista que no contiene solo cadenas
        with self.assertRaises(ValueError):
            handler.download_fastas([123, 'PDB2'])

    @patch('os.path.exists')
    @patch('concurrent.futures.ThreadPoolExecutor')
    @patch('protein_data_handler.fasta.FastaHandler.download_and_store_fasta')
    def test_download_fastas_with_existing_files(self, mock_download_and_store_fasta, mock_executor, mock_exists):
        # Configuración para simular la existencia de un archivo para el primer pdb_id
        def side_effect(path):
            return 'PDB1.fasta' in path

        mock_exists.side_effect = side_effect
        mock_download_and_store_fasta.return_value = ('path', [('PDB2', 'chain1', 'id1', 'sequence1')])
        session_mock = MagicMock()
        handler = FastaHandler(session_mock, './data', './output')

        pdb_ids = ['PDB1', 'PDB2']
        handler.download_fastas(pdb_ids)

    @patch('os.path.exists')
    @patch('os.makedirs')
    @patch('requests.get')
    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    def test_download_and_store_fasta(self, mock_open, mock_get, mock_makedirs, mock_exists):
        mock_exists.return_value = False
        mock_get.return_value = MagicMock(status_code=200, text='fasta data')
        session_mock = MagicMock()
        handler = FastaHandler(session_mock, './data', './output')
        file_path, chains = handler.download_and_store_fasta('PDB1')
        mock_open.assert_called()
        self.assertIsNotNone(file_path)
        self.assertIsNotNone(chains)

        # Test para manejar RequestException
        mock_get.side_effect = RequestException()
        file_path, chains = handler.download_and_store_fasta('PDB1')
        self.assertIsNone(file_path)
        self.assertIsNone(chains)

        # Test para manejar IOError
        mock_open.side_effect = IOError()
        file_path, chains = handler.download_and_store_fasta('PDB1')
        self.assertIsNone(file_path)
        self.assertIsNone(chains)

    @patch('os.path.exists')
    @patch('os.makedirs')
    @patch('requests.get')
    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    def test_download_and_store_fasta_invalid_input(self, mock_open, mock_get, mock_makedirs, mock_exists):
        mock_exists.return_value = False
        session_mock = MagicMock()
        handler = FastaHandler(session_mock, './data', './output')

        # Prueba con un argumento no válido (no es una cadena)
        with self.assertRaises(ValueError):
            handler.download_and_store_fasta(123)  # Uso de un número en lugar de una cadena

        # Otra prueba podría ser pasar un objeto o una lista
        with self.assertRaises(ValueError):
            handler.download_and_store_fasta(['PDB1'])  # Uso de una lista en lugar de una cadena

    @patch('os.path.exists')
    @patch('requests.get')
    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    def test_download_and_store_fasta_io_error(self, mock_open, mock_get, mock_exists):
        mock_exists.return_value = False
        mock_get.return_value = MagicMock(status_code=200, text='fasta data')

        # Configura mock_open para lanzar IOError al intentar escribir en un archivo
        mock_open.side_effect = IOError("IO error")

        session_mock = MagicMock()
        handler = FastaHandler(session_mock, './data', './output')

        file_path, chains = handler.download_and_store_fasta('PDB1')

        # Verificar que se capturó y se registró el IOError
        self.assertIsNone(file_path)
        self.assertIsNone(chains)

    @patch('protein_data_handler.fasta.cd_hit')
    @patch('protein_data_handler.fasta.read_clstr')
    def test_cluster_fastas(self, mock_read_clstr, mock_cd_hit):
        # DataFrame simulado que podría ser retornado por read_clstr
        simulated_df = pd.DataFrame({
            'cluster': ['cluster1'],
            'identifier': ['PDB1_chain1'],
            'is_representative': [1],
            'size': [100],
            'identity': [0.99]
        })
        mock_read_clstr.return_value = simulated_df

        # Simulamos la función cd_hit para que no ejecute un proceso externo
        mock_cd_hit.return_value = None

        session_mock = MagicMock()
        handler = FastaHandler(session_mock, './data', './output')

        # Ejecutar el método cluster_fastas
        df_result = handler.cluster_fastas('input_file')
        mock_cd_hit.assert_called()
        mock_read_clstr.assert_called()
        self.assertIsNotNone(df_result)

    @patch('os.path.join')
    @patch('os.path.isfile')
    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    def test_merge_fastas(self, mock_open, mock_isfile, mock_join):
        session_mock = MagicMock()
        handler = FastaHandler(session_mock, './data', './output')
        mock_isfile.return_value = True
        pdb_ids = ['PDB1', 'PDB2']
        handler.merge_fastas(pdb_ids, 'merged')
        mock_open.assert_called()

    @patch('os.path.join')
    @patch('os.path.isfile')
    @patch('builtins.open', new_callable=unittest.mock.mock_open)
    @patch('logging.warning')
    def test_merge_fastas_file_not_found(self, mock_logging_warning, mock_open, mock_isfile, mock_join):
        session_mock = MagicMock()
        handler = FastaHandler(session_mock, './data', './output')

        # Configura mock_isfile para simular que uno de los archivos no existe
        def side_effect(file_path):
            return 'PDB2.fasta' in file_path

        mock_isfile.side_effect = side_effect

        pdb_ids = ['PDB1', 'PDB2']
        merge_name = 'merged'
        handler.merge_fastas(pdb_ids, merge_name)

