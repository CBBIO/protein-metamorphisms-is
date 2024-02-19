from unittest.mock import patch, MagicMock, call
import unittest

# Asegúrate de que la importación de PDBExtractor sea correcta según tu estructura de paquetes
from protein_metamorphisms_is.information_system.pdb import PDBExtractor

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
        with patch('protein_metamorphisms_is.information_system.base.extractor.setup_logger'), \
                patch('protein_metamorphisms_is.information_system.base.extractor.create_engine'), \
                patch('protein_metamorphisms_is.information_system.base.extractor.sessionmaker', return_value=MagicMock()) as mock_sessionmaker:
            self.mock_session = mock_sessionmaker.return_value
            self.extractor = PDBExtractor(self.conf)


    def test_load_pdb_ids(self):
        # Configura el mock de la sesión para simular el comportamiento de la consulta y su resultado
        mock_query = self.mock_session.query.return_value
        mock_filter = mock_query.filter.return_value
        mock_filter.all.return_value = [MagicMock(pdb_id='1XYZ', resolution=1.5), MagicMock(pdb_id='2ABC', resolution=2.0)]

        with patch('sqlalchemy.sql.operators.or_'):

            pdb_references = self.extractor.load_pdb_ids()

            # Verifica que se devuelvan las referencias correctas
            self.assertEqual(len(pdb_references), 2)
            self.assertEqual(pdb_references[0].pdb_id, '1XYZ')
            self.assertEqual(pdb_references[1].pdb_id, '2ABC')
