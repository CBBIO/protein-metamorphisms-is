from protein_data_handler.uniprot import (
    descargar_registro,
    cargar_codigos_acceso,
    extraer_entradas,
    almacenar_entrada,
)
import unittest
from unittest.mock import patch, MagicMock
import requests
from concurrent.futures import ThreadPoolExecutor, Future


class TestCargarCodigosAcceso(unittest.TestCase):
    @patch("requests.get")
    def test_cargar_codigos_acceso_exitoso(self, mock_get):
        # Configura el mock de requests.get para simular una respuesta exitosa
        mock_get.return_value.ok = True
        mock_get.return_value.text = "ABC123\nDEF456\nGHI789"
        mock_session = MagicMock()

        # Llama a la función con los mocks
        cargar_codigos_acceso("criterio_de_busqueda", 3, mock_session)

        # Verifica que se realizaron las llamadas esperadas a la sesión de la base de datos
        self.assertEqual(mock_session.add.call_count, 3)
        mock_session.commit.assert_called_once()

    @patch("requests.get")
    def test_cargar_codigos_acceso_falla_api(self, mock_get):
        # Configura el mock de requests.get para simular un error
        mock_get.return_value.raise_for_status.side_effect = (
            requests.exceptions.HTTPError
        )

        mock_session = MagicMock()

        # Llama a la función esperando un error
        cargar_codigos_acceso("criterio_de_busqueda", 3, mock_session)

        # Verifica que se realizó un rollback en la sesión
        mock_session.rollback.assert_called_once()


class TestDescargarRegistro(unittest.TestCase):
    @patch("Bio.ExPASy.get_sprot_raw")
    @patch("Bio.SwissProt.read")
    def test_successful_download(self, mock_read, mock_get_sprot_raw):
        mock_handle = MagicMock()
        mock_record = MagicMock()
        mock_get_sprot_raw.return_value = mock_handle
        mock_read.return_value = mock_record

        result = descargar_registro("valid_accession_code")
        self.assertEqual(result, mock_record)

    @patch("Bio.ExPASy.get_sprot_raw")
    def test_invalid_accession_code(self, mock_get_sprot_raw):
        mock_get_sprot_raw.side_effect = Exception("Invalid accession code")

        result = descargar_registro("invalid_code")
        self.assertIsNone(result)

    @patch("Bio.ExPASy.get_sprot_raw")
    def test_exception_handling(self, mock_get_sprot_raw):
        mock_get_sprot_raw.side_effect = Exception("Network error")

        result = descargar_registro("some_code")
        self.assertIsNone(result)


class TestExtraerEntradas(unittest.TestCase):
    @patch("protein_data_handler.uniprot.descargar_registro")
    @patch("protein_data_handler.uniprot.almacenar_entrada")
    def test_extraer_entradas_exitoso(self, mock_almacenar, mock_descargar):
        mock_session = MagicMock()
        mock_session.query.return_value.all.return_value = [
            MagicMock(accession_code="ABC123"),
            MagicMock(accession_code="DEF456"),
        ]

        # Simular que descargar_registro devuelve un resultado
        mock_descargar.return_value = "datos_de_prueba"

        # Llamar a la función
        extraer_entradas(mock_session)

        # Verificar que se llamó a almacenar_entrada para cada descarga exitosa
        self.assertEqual(mock_almacenar.call_count, 2)

    @patch("protein_data_handler.uniprot.descargar_registro")
    @patch("protein_data_handler.uniprot.almacenar_entrada")
    def test_extraer_entradas_con_error(self, mock_almacenar, mock_descargar):
        mock_session = MagicMock()
        mock_session.query.return_value.all.return_value = [
            MagicMock(accession_code="ABC123")
        ]
        mock_almacenar.side_effect = Exception(
            "Fallo al insertar entrada error"
        )
        # Llamar a la función
        extraer_entradas(mock_session)

        # Como hubo un error, almacenar_entrada no debería ser llamado
        mock_session.add.assert_not_called()


class TestAlmacenarEntrada(unittest.TestCase):
    def test_creacion_nueva_proteina(self):
        mock_data = MagicMock()
        mock_data.entry_name = "nueva_proteina"
        mock_data.accessions = ["accession_1", "accession_2"]

        mock_session = MagicMock()
        mock_session.query.return_value.filter_by.return_value.first.return_value = (
            None
        )

        almacenar_entrada(mock_data, mock_session)

        # Verificar que se añadió la nueva proteína
        self.assertTrue(mock_session.add.called)

    def test_procesamiento_referencias_pdb(self):
        mock_data = MagicMock()
        mock_data.entry_name = "proteina_test"
        mock_data.cross_references = [
            ("PDB", "1234", "método_test", "resolución_test", "A=10-444, D=59-430")
        ]
        mock_session = MagicMock()
        mock_session.query.return_value.filter_by.return_value.first.return_value = (
            None
        )

        almacenar_entrada(mock_data, mock_session)

        # Verificar que se añadió la referencia PDB
        self.assertTrue(mock_session.add.called)

    def test_procesamiento_terminos_go(self):
        mock_data = MagicMock()
        mock_data.entry_name = "proteina_test"
        mock_data.cross_references = [
            ("GO", "GO:0001234", "categoría_test:descripción_test")
        ]

        mock_session = MagicMock()
        mock_session.query.return_value.filter_by.return_value.first.return_value = (
            None
        )

        almacenar_entrada(mock_data, mock_session)

        # Verificar que se añadió el término GO
        self.assertTrue(mock_session.add.called)

    def test_actualizacion_proteina_existente(self):
        mock_data = MagicMock()
        mock_data.entry_name = "proteina_existente"
        mock_session = MagicMock()
        mock_session.query.return_value.filter_by.return_value.first.return_value = (
            MagicMock()
        )

        almacenar_entrada(mock_data, mock_session)

        # Verificar que la proteína existente se actualizó
        self.assertTrue(mock_session.add.called)

    def test_manejo_excepciones(self):
        mock_session = MagicMock()
        mock_session.query.side_effect = Exception("Error de prueba")

        almacenar_entrada(MagicMock(), mock_session)

        # Verificar que se hizo rollback
        mock_session.rollback.assert_called_once()

    @patch("protein_data_handler.models.uniprot.PDBReference")
    @patch("protein_data_handler.models.uniprot.Chain")
    def test_procesamiento_cadenas_pdb(self, mock_chain_class, mock_pdb_reference_class):
        # Crear un mock para los datos de entrada y la sesión
        mock_data = MagicMock()
        mock_data.entry_name = "proteina_test"
        mock_data.cross_references = [
            ("PDB", "1234", "método_test", "resolución_test", "A=10-444, D=59-430")
        ]

        mock_session = MagicMock()

        # Configurar el comportamiento esperado de la sesión para PDBReference
        mock_pdb_reference_instance = MagicMock()
        mock_pdb_reference_instance.id = 1


        mock_chain_instance = MagicMock()
        mock_chain_class.return_value = mock_chain_instance

        mock_session.query.return_value.filter_by.return_value.first.side_effect = [
            mock_pdb_reference_instance,  # Para PDBReference
            None, 
            mock_chain_instance,
            None,

        ]



        # Llamar a la función almacenar_entrada
        almacenar_entrada(mock_data, mock_session)


        # self.assertEqual(mock_session.query.return_value.filter_by.return_value.first.call_count,2)


        # Verificar que se buscó la referencia PDB y se procesaron las cadenas
        expected_pdb_calls = [((mock_pdb_reference_class, 'pdb_id'), {'pdb_id': '1234'})]
        expected_chain_calls = [
            ((mock_chain_class, 'pdb_reference_id', 'chain'), {'pdb_reference_id': 1, 'chain': 'A'}),
            ((mock_chain_class, 'pdb_reference_id', 'chain'), {'pdb_reference_id': 1, 'chain': 'D'}),
        ]
        # self.assertEqual(mock_session.query.call_args_list, expected_pdb_calls + expected_chain_calls)

        # Verificar que se añadieron las cadenas a la sesión si no existían
        # self.assertEqual(mock_session.add.call_count, 2)


if __name__ == "__main__":
    unittest.main()
