import os
import unittest
from unittest.mock import patch, MagicMock, call

import requests
from requests import RequestException

from protein_data_handler.fasta import FastaHandler


class TestFastaHandler(unittest.TestCase):

    def setUp(self):
        self.session_mock = MagicMock()  # Simula la sesión de SQLAlchemy
        self.data_dir = './tests/data/FASTA'
        self.output_dir = './tests/data/FASTA/output'
        self.downloader = FastaHandler(self.session_mock, self.data_dir, self.output_dir)

    def test_init(self):
        self.assertEqual(self.downloader.session, self.session_mock)

    @patch('os.makedirs')
    @patch('os.path.exists')
    def test_creacion_directorio_si_no_existe(self, mock_exists, mock_makedirs):
        # Configurar el mock para simular que el directorio no existe
        mock_exists.return_value = False

        # Inicializar FastaHandler
        FastaHandler('sesion_falsa', './tests/data/FASTA', './tests/data/FASTA/output')

        # Verificar que os.makedirs fue llamado
        expected_calls = [call('./tests/data/FASTA', exist_ok=True),
                          call('./tests/data/FASTA/output', exist_ok=True)]
        mock_makedirs.assert_has_calls(expected_calls)

    def test_download_fastas_invalid_input(self):
        with self.assertRaises(ValueError):
            self.downloader.download_fastas("no-list")

    @patch('protein_data_handler.fasta.requests.get')
    def test_download_fasta_successful(self, mock_get):
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "fake-fasta-content"
        mock_get.return_value = mock_response
        with patch('builtins.open', unittest.mock.mock_open()) as mock_file:
            self.downloader.download_fasta("PDBID")
            mock_file.assert_called_with("./tests/data/FASTA/PDBID.fasta", "w")
            mock_file().write.assert_called_with("fake-fasta-content")

    @patch('protein_data_handler.fasta.requests.get')
    @patch('protein_data_handler.fasta.logging.error')
    def test_download_fasta_http_error(self, mock_logging_error, mock_get):
        mock_get.side_effect = requests.exceptions.RequestException("HTTP error")

        self.downloader.download_fasta("PDBID")
        mock_logging_error.assert_called_with("Error al descargar FASTA para PDBID: HTTP error")

    @patch('protein_data_handler.fasta.requests.get')
    @patch('protein_data_handler.fasta.logging.error')
    def test_download_fasta_io_error(self, mock_logging_error, mock_get):
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "fake-fasta-content"
        mock_get.return_value = mock_response

        with patch('builtins.open', unittest.mock.mock_open()) as mock_file:
            mock_file.side_effect = IOError("IO error")
            self.downloader.download_fasta("PDBID")
            mock_logging_error.assert_called_with("Error al escribir el archivo para PDBID: IO error")

    @patch('protein_data_handler.fasta.FastaHandler.download_fasta')
    def test_download_fastas(self, mock_download_fasta):
        pdb_ids = ["PDB1", "PDB2", "PDB3"]
        self.downloader.download_fastas(pdb_ids)
        self.assertEqual(mock_download_fasta.call_count, len(pdb_ids))

    def test_download_fasta_invalid_pdb_id(self):
        with self.assertRaises(ValueError):
            self.downloader.download_fasta(123)

    @patch('os.path.exists')
    @patch('os.makedirs')
    def test_download_fasta_directory_creation(self, mock_makedirs, mock_exists):
        mock_exists.return_value = None
        mock_makedirs.return_value = None

        self.downloader.download_fasta("PDBID")
        mock_exists.return_value = True
        mock_makedirs.return_value = True

        self.downloader.download_fasta("PDBID")

    @patch('os.path.isfile')
    @patch('builtins.open', new_callable=unittest.mock.mock_open, read_data="contenido_fasta")
    @patch('protein_data_handler.fasta.FastaHandler.download_fastas')
    def test_merge_fastas(self, mock_download_fastas, mock_open, mock_isfile):
        # Configura los mocks
        mock_isfile.side_effect = lambda x: x.endswith('.fasta')
        pdb_ids = ['PDB1', 'PDB2', 'PDB3']
        merge_name = 'merged'

        # Ejecuta la función
        self.downloader.merge_fastas(pdb_ids, merge_name)

        # Verifica que se abran los archivos correctos y se escriba en el archivo de salida
        expected_file_calls = [call(os.path.join(self.downloader.data_dir, f'{pdb_id}.fasta'), 'r') for pdb_id in
                               pdb_ids]
        expected_file_calls.append(call(os.path.join(self.downloader.output_dir, f'{merge_name}.fasta'), 'w'))
        mock_open.assert_has_calls(expected_file_calls, any_order=True)

    @patch('os.path.isfile')
    @patch('builtins.open', new_callable=unittest.mock.mock_open, read_data="contenido_fasta")
    @patch('protein_data_handler.fasta.FastaHandler.download_fastas')
    def test_merge_fastas_fichero_inexistente(self, mock_download_fastas, mock_open, mock_isfile):
        # Configura los mocks
        mock_isfile.side_effect = lambda x: x.endswith('.fasta')
        pdb_ids = ['PDB1', 'PDB2', 'PDB3']
        merge_name = 'merged'

        # Ejecuta la función
        self.downloader.merge_fastas(pdb_ids, merge_name)

        # Verifica que se abran los archivos correctos y se escriba en el archivo de salida
        expected_file_calls = [call(os.path.join(self.downloader.data_dir, f'{pdb_id}.fasta'), 'r') for pdb_id in
                               pdb_ids]
        expected_file_calls.append(call(os.path.join(self.downloader.output_dir, f'{merge_name}.fasta'), 'w'))
        mock_open.assert_has_calls(expected_file_calls, any_order=True)

    @patch('os.path.isfile')
    @patch('builtins.open', new_callable=unittest.mock.mock_open, read_data="contenido_fasta")
    @patch('protein_data_handler.fasta.FastaHandler.download_fastas')
    @patch('logging.info')
    def test_merge_fastas_with_missing_files(self, mock_logging_info, mock_download_fastas, mock_open, mock_isfile):
        # Configura los mocks para simular que algunos archivos no existen
        mock_isfile.side_effect = lambda filepath: 'PDB2.fasta' in filepath
        pdb_ids = ['PDB1', 'PDB2', 'PDB3']
        merge_name = 'merged'

        # Ejecuta la función
        self.downloader.merge_fastas(pdb_ids, merge_name)

        # Verifica que se llame a download_fastas para los archivos faltantes
        mock_download_fastas.assert_called_with(['PDB1', 'PDB3'])

        # Verifica que se registre la información de descarga
        mock_logging_info.assert_called_with("Descargando archivos FASTA faltantes.")

        # Verifica que se abran los archivos correctos y se escriba en el archivo de salida
        expected_file_calls = [call(os.path.join(self.data_dir, 'PDB2.fasta'), 'r')]
        expected_file_calls.append(call(os.path.join(self.output_dir, f'{merge_name}.fasta'), 'w'))
        mock_open.assert_has_calls(expected_file_calls, any_order=True)

    @patch('os.path.isfile')
    @patch('logging.warning')
    def test_merge_fastas_file_not_found(self, mock_logging_warning, mock_isfile):
        # Configura los mocks para simular que ningún archivo existe
        mock_isfile.return_value = False
        pdb_ids = ['PDB1']
        merge_name = 'merged'

        # Ejecuta la función
        self.downloader.merge_fastas(pdb_ids, merge_name)

        # Verifica que se registre una advertencia para el archivo no encontrado
        mock_logging_warning.assert_called_with(f"Archivo no encontrado: ./tests/data/FASTA/PDB1.fasta")


if __name__ == '__main__':
    unittest.main()
