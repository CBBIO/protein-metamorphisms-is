import unittest
from unittest.mock import patch

from protein_data_handler.pdb import download_entire_pdb


class TestDownloadEntirePDB(unittest.TestCase):

    @patch('protein_data_handler.pdb.PDBList')
    def test_download_entire_pdb(self, mock_pdblist):
        # Configura el mock
        mock_instance = mock_pdblist.return_value
        mock_instance.download_entire_pdb.return_value = None

        # Llama a la función con parámetros de prueba
        server = "https://testserver.org/"
        pdb = "./test_data"
        file_format = 'mmCif'
        download_entire_pdb(server, pdb, file_format)

        # Asegúrate de que se llamó a PDBList con los parámetros correctos
        mock_pdblist.assert_called_with(server=server, pdb=pdb)

        # Verifica que se llamó a download_entire_pdb con el formato de archivo correcto
        mock_instance.download_entire_pdb.assert_called_with(file_format=file_format)

