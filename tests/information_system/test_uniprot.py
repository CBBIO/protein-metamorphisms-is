import unittest
from concurrent.futures import Future
from http.client import HTTPException
from unittest.mock import patch, MagicMock, call

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
        with patch('protein_metamorphisms_is.information_system.base.extractor.setup_logger') as mock_logger, \
                patch('protein_metamorphisms_is.sql.base.database_manager.create_engine'), \
                patch('protein_metamorphisms_is.sql.base.database_manager.sessionmaker',
                      return_value=MagicMock()):
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
    def test_load_access_codes_with_new_accessions(self, mock_get):
        # Configura el mock de requests.get para simular una respuesta exitosa
        mock_response = MagicMock(status_code=200, text="P12345\nP67890")
        mock_get.return_value = mock_response

        # Configura el mock de session.query().scalar() para simular que los códigos de acceso no existen
        self.extractor.session.query.return_value.scalar.side_effect = [False, False]

        # Ejecuta la función bajo prueba
        self.extractor.load_access_codes("example search criteria", 2)

        # Verifica que requests.get fue llamado correctamente
        mock_get.assert_called_once()

        # Verifica que session.add fue llamado para cada nuevo código de acceso
        self.assertEqual(self.extractor.session.add.call_count, 2)

        # Verifica que session.commit fue llamado para confirmar los cambios
        self.extractor.session.commit.assert_called_once()

    @patch('requests.get')
    def test_load_access_codes_exception_handling(self, mock_get):
        # Configura el mock de requests.get para lanzar una excepción
        mock_get.side_effect = requests.exceptions.HTTPError("Error de solicitud")

        # Ejecuta la función bajo prueba dentro de un bloque try-except para capturar la excepción
        try:
            self.extractor.load_access_codes("example search criteria", 2)
        except requests.exceptions.HTTPError:
            pass  # Aquí simplemente pasamos ya que estamos probando el manejo de la excepción dentro de la función

        # Verifica que se llamó al rollback de la sesión debido a la excepción
        self.extractor.session.rollback.assert_called_once()

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

    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor.download_record')
    @patch('protein_metamorphisms_is.information_system.uniprot.UniProtExtractor.store_entry')
    def test_extract_entries(self, mock_store_entry, mock_download_record):
        # Configura mock_download_record para devolver datos, None, y luego lanzar una excepción
        mock_download_record.side_effect = [MagicMock(), None, Exception("Failed to download")]

        # Configura la respuesta de la base de datos simulada con tres accesiones
        self.extractor.session.query.return_value.all.return_value = [
            MagicMock(accession_code="ABC123"),
            MagicMock(accession_code="DEF456"),
            MagicMock(accession_code="GHI789"),
        ]

        self.extractor.extract_entries()

        # Verifica que download_record fue llamado tres veces, una para cada accession_code
        expected_calls = [call('ABC123'), call('DEF456'), call('GHI789')]
        mock_download_record.assert_has_calls(expected_calls, any_order=True)

        # Verifica que store_entry fue llamado solo una vez
        mock_store_entry.assert_called_once()

        # Verifica que el logger registró el error para la tercera llamada que lanza una excepción
        self.extractor.logger.error.assert_called()

    def test_store_entry(self):
        # Crea un registro SwissProt.Record simulado
        mock_record = MagicMock()
        mock_record.entry_name = 'SampleEntry'
        mock_record.data_class = 'SampleClass'
        # Configura los demás campos necesarios de mock_record según tu implementación

        # Ejecuta store_entry
        self.extractor.store_entry(mock_record)

        # Verifica que Session.add fue llamado
        self.extractor.session.add.assert_called()
        self.extractor.session.commit.assert_called()

    def test_store_entry_with_new_accession_and_pdb_reference(self):
        # Prepara un registro SwissProt.Record simulado con los datos necesarios
        mock_record = MagicMock()
        mock_record.entry_name = 'SampleEntry'
        mock_record.data_class = 'SampleClass'
        mock_record.accessions = ['A12345']
        mock_record.cross_references = [('PDB', 'P12345', 'Method', '2.0')]
        # Añade más configuraciones al mock_record según sea necesario

        # Configura los mocks para simular la verificación de existencia y la creación de nuevas entradas
        self.extractor.session.scalar.side_effect = [False, False]

        # Ejecuta store_entry
        self.extractor.store_entry(mock_record)

        # Verifica que se agregaron las entradas y referencias nuevas a la sesión
        self.assertEqual(self.extractor.session.add.call_count,
                         1)  # Cambia este número según el número esperado de llamadas a add()
        # Verifica que se llamó a commit para confirmar los cambios
        self.extractor.session.commit.assert_called_once()

        # Verifica que se haya llamado correctamente a la consulta de existencia
        self.assertEqual(self.extractor.session.query.call_count,
                         4)  # Ajusta este número según el número esperado de llamadas a query()

    def test_store_entry_updates_existing_protein(self):
        mock_record = MagicMock()
        mock_record.entry_name = 'ExistingEntry'
        mock_record.data_class = 'UpdatedClass'
        # Configura los demás campos necesarios de mock_record según tu implementación

        # Simula que la entrada de proteína ya existe
        self.extractor.session.query.return_value.scalar.return_value = True
        self.extractor.session.query.return_value.filter_by.return_value.first.return_value = MagicMock()

        self.extractor.store_entry(mock_record)

        # Verifica que se actualizó la entrada existente sin intentar crear una nueva
        self.extractor.session.add.assert_called()
        self.extractor.session.commit.assert_called()

    def test_store_entry_exception_handling(self):
        mock_record = MagicMock()
        mock_record.entry_name = 'FaultyEntry'
        # Configura los demás campos necesarios de mock_record según tu implementación

        # Simula una excepción al intentar hacer commit
        self.extractor.session.commit.side_effect = HTTPException("Commit failed")

        with patch.object(self.extractor.logger, 'error') as mock_log_error:
            self.extractor.store_entry(mock_record)

            # Verifica que se llamó a rollback debido a la excepción
            self.extractor.session.rollback.assert_called_once()

            # Verifica que el logger registró el error
            self.extractor.logger.error.assert_called_with("Error while dumping the entry: Commit failed")

            # Asegura que la sesión se cierra correctamente
            self.extractor.session.close.assert_called_once()

    def test_store_entry_creates_new_protein_when_not_exists(self):
        # Crea un registro SwissProt.Record simulado con los datos necesarios
        mock_record = MagicMock()
        mock_record.entry_name = 'NewEntry'
        mock_record.data_class = 'SampleClass'
        mock_record.molecule_type = 'Protein'
        mock_record.sequence_length = 123
        mock_record.created = ['2022-01-01']
        mock_record.sequence = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF...'
        mock_record.sequence_update = ['2022-01-02']
        mock_record.annotation_update = ['2022-01-03']
        mock_record.description = 'A sample protein'
        mock_record.gene_name = 'GENE1'
        mock_record.organism = 'Sample Organism'
        mock_record.organelle = 'Cell'
        mock_record.organism_classification = ['Eukaryota', 'Metazoa']
        mock_record.taxonomy_id = ['9606']
        mock_record.host_organism = ['Host Organism']
        mock_record.host_taxonomy_id = ['9607']
        mock_record.comments = ['Function: Sample function.']
        mock_record.keywords = ['Keyword1', 'Keyword2']
        mock_record.protein_existence = 'Evidence at protein level'
        mock_record.seqinfo = (123, 'A123B456', 'V1')
        mock_record.accessions = ['A12345']
        mock_record.cross_references = [('PDB', 'P12345', 'Method', '2.0')]

        self.extractor.session.query.return_value.scalar.side_effect = [False,
                                                                        False, False]

        # Ejecuta store_entry
        self.extractor.store_entry(mock_record)

        # Verifica que se creó una nueva entrada de proteína y se agregó a la sesión
        self.extractor.session.add.assert_called()
        # Verifica que se llamó a commit para confirmar los cambios
        self.extractor.session.commit.assert_called_once()


if __name__ == '__main__':
    unittest.main()
