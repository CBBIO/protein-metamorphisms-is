import unittest
from concurrent.futures import Future
from unittest.mock import patch, MagicMock

import requests

from protein_metamorphisms_is.information_system.uniprot import UniProtExtractor

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
            "max_workers": 2  # Asegúrate de definir esto si tu clase lo requiere
        }
        self.extractor = UniProtExtractor(self.conf)

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor.extract_entries')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor.load_access_codes')
    def test_start(self, mock_load_access_codes, mock_extract_entries):
        # Caso de prueba para un flujo normal sin excepciones
        self.extractor.start()
        mock_load_access_codes.assert_called_once_with("(structure_3d:true)", 100)
        mock_extract_entries.assert_called_once()

        # Resetea los mocks para el siguiente caso de prueba
        mock_load_access_codes.reset_mock()
        mock_extract_entries.reset_mock()

        # Caso de prueba para manejo de excepciones
        mock_load_access_codes.side_effect = Exception("Error de prueba")
        with patch.object(self.extractor.logger, 'error') as mock_log_error:
            self.extractor.start()
            mock_log_error.assert_called_with("Error during extraction process: Error de prueba")

        # Asegúrate de que, a pesar del error, se intentó comenzar el proceso de extracción
        mock_load_access_codes.assert_called_once()
        # `extract_entries` no debe ser llamado debido al error en `load_access_codes`
        mock_extract_entries.assert_not_called()

    @patch('requests.get')
    def test_load_access_codes(self, mock_get):
        # Simula una respuesta exitosa de requests.get
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "P12345\nP67890"
        mock_get.return_value = mock_response

        # Simula las interacciones con la base de datos aquí
        with patch('sqlalchemy.orm.Session.query') as mock_query:
            # Agrega más lógica de simulación según sea necesario para tu caso de uso

            self.extractor.load_access_codes("example search criteria", 2)
            # Asegúrate de que requests.get fue llamado correctamente
            mock_get.assert_called_once()
            # Añade más aserciones para verificar el comportamiento de tu función

    @patch('requests.get')
    @patch('sqlalchemy.orm.Session.add')
    @patch('sqlalchemy.orm.Session.commit')
    @patch('sqlalchemy.orm.session.Session.query')
    def test_load_access_codes_with_new_accessions(self, mock_query, mock_commit, mock_add, mock_get):
        # Configura el mock de requests.get para simular una respuesta exitosa
        mock_response = MagicMock(status_code=200, text="P12345\nP67890")
        mock_get.return_value = mock_response

        # Configura el mock de session.query().scalar() para simular que los códigos de acceso no existen
        mock_query.return_value.scalar.return_value = False

        # Ejecuta la función bajo prueba
        self.extractor.load_access_codes("example search criteria", 2)

        # Verifica que requests.get fue llamado correctamente
        mock_get.assert_called_once()

        # Verifica que session.add fue llamado para cada nuevo código de acceso
        self.assertEqual(mock_add.call_count, 2)

        # Verifica que session.commit fue llamado para confirmar los cambios
        mock_commit.assert_called_once()

    @patch('requests.get')
    @patch('sqlalchemy.orm.session.Session.rollback')
    @patch('sqlalchemy.orm.session.Session.query')
    def test_load_access_codes_exception_handling(self, mock_query, mock_rollback, mock_get):
        # Configura el mock de requests.get para lanzar una excepción
        mock_get.side_effect = requests.exceptions.HTTPError("Error de solicitud")

        # Ejecuta la función bajo prueba dentro de un bloque try-except para capturar la excepción
        try:
            self.extractor.load_access_codes("example search criteria", 2)
        except requests.exceptions.HTTPError:
            pass  # Aquí simplemente pasamos ya que estamos probando el manejo de la excepción dentro de la función

        # Verifica que se llamó al rollback de la sesión debido a la excepción
        mock_rollback.assert_called_once()

        # Verifica que el logger registró el error
        with patch.object(self.extractor.logger, 'error') as mock_log_error:
            self.extractor.load_access_codes("example search criteria", 2)
            mock_log_error.assert_called_with("Error: Error de solicitud")
    @patch('Bio.ExPASy.get_sprot_raw')
    @patch('Bio.SwissProt.read')
    def test_download_record(self, mock_read, mock_get_sprot_raw):
        # Simula una respuesta exitosa de ExPASy.get_sprot_raw y SwissProt.read
        mock_get_sprot_raw.return_value = MagicMock()
        mock_record = MagicMock()
        mock_read.return_value = mock_record

        result = self.extractor.download_record("P12345")
        self.assertIsNotNone(result)
        # Verifica que el resultado es el esperado
        mock_get_sprot_raw.assert_called_once_with("P12345")
        mock_read.assert_called_once()

    @patch('Bio.ExPASy.get_sprot_raw')
    def test_download_record_exception_handling(self, mock_get_sprot_raw):
        # Configura el mock para lanzar una excepción
        mock_get_sprot_raw.side_effect = Exception("Error de prueba")

        # Configura un mock para el logger
        with patch.object(self.extractor.logger, 'error') as mock_log_error:
            # Ejecuta el método bajo prueba
            result = self.extractor.download_record("P12345")

            # Verifica que el resultado sea None debido a la excepción
            self.assertIsNone(result)

            # Verifica que el logger haya registrado el error
            mock_log_error.assert_called_with("Error downloading the entry P12345: Error de prueba")

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor.download_record',
           return_value=MagicMock())
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor.store_entry')
    @patch('sqlalchemy.orm.Session.query')
    def test_extract_entries(self, mock_query, mock_store_entry, mock_download_record):
        # Primero, mockeamos el retorno del método .all() de la consulta
        mock_accession = MagicMock(accession_code='P12345')
        mock_query.return_value.all.return_value = [mock_accession]

        # Ahora ejecutamos extract_entries
        self.extractor.extract_entries()

        # Verificaciones
        mock_download_record.assert_called_once_with('P12345')
        mock_store_entry.assert_called()

    from unittest.mock import patch, MagicMock
    from concurrent.futures import Future
    import logging

    @patch('sqlalchemy.orm.session.Session.query')
    def test_extract_entries_exception_handling(self, mock_query):
        # Configura la respuesta de la base de datos simulada
        mock_accession = MagicMock(accession_code='P12345')
        mock_query.return_value.all.return_value = [mock_accession]

        # Configura el mock del ThreadPoolExecutor y future.result para lanzar una excepción
        mock_future = MagicMock(spec=Future)
        mock_future.result.side_effect = Exception("Error durante la descarga")

        # Configura un mock para el logger
        with patch.object(self.extractor.logger, 'error') as mock_log_error:
            self.extractor.extract_entries()

            mock_log_error.assert_called()


    @patch('sqlalchemy.orm.Session.add')
    @patch('sqlalchemy.orm.Session.commit')
    def test_store_entry(self, mock_commit, mock_add):
        # Crea un registro SwissProt.Record simulado
        mock_record = MagicMock()
        mock_record.entry_name = 'SampleEntry'
        mock_record.data_class = 'SampleClass'
        # Configura los demás campos necesarios de mock_record según tu implementación

        # Ejecuta store_entry
        self.extractor.store_entry(mock_record)

        # Verifica que Session.add fue llamado
        mock_add.assert_called()
        # Verifica que Session.commit fue llamado para confirmar los cambios
        mock_commit.assert_called()

    @patch('sqlalchemy.orm.session.Session.add')
    @patch('sqlalchemy.orm.session.Session.commit')
    @patch('sqlalchemy.orm.session.Session.query')
    def test_store_entry_with_new_accession_and_pdb_reference(self, mock_query, mock_commit, mock_add):
        # Prepara un registro SwissProt.Record simulado con los datos necesarios
        mock_record = MagicMock()
        mock_record.entry_name = 'SampleEntry'
        mock_record.data_class = 'SampleClass'
        mock_record.accessions = ['A12345']
        mock_record.cross_references = [('PDB', 'P12345', 'Method', '2.0')]
        # Añade más configuraciones al mock_record según sea necesario

        # Configura los mocks para simular la verificación de existencia y la creación de nuevas entradas
        mock_query.scalar.side_effect = [False, False]

        # Ejecuta store_entry
        self.extractor.store_entry(mock_record)

        # Verifica que se agregaron las entradas y referencias nuevas a la sesión
        self.assertEqual(mock_add.call_count, 1)  # Cambia este número según el número esperado de llamadas a add()
        # Verifica que se llamó a commit para confirmar los cambios
        mock_commit.assert_called_once()

        # Verifica que se haya llamado correctamente a la consulta de existencia
        self.assertEqual(mock_query.call_count, 4)  # Ajusta este número según el número esperado de llamadas a query()


if __name__ == '__main__':
    unittest.main()
