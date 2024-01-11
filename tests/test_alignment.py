import unittest
from unittest.mock import patch, MagicMock
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from protein_data_handler.alignment import UniProtPDBMapping
from protein_data_handler.sql.model import PDBReference

class TestUniProtPDBMapping(unittest.TestCase):

    @patch('Bio.AlignIO.read')
    @patch('subprocess.Popen')
    def test_alinear_secuencias_mafft(self, mock_popen, mock_alignio_read):
        # Configurar mocks
        process_mock = MagicMock()
        attrs = {'communicate.return_value': (">Seq1\nAAA\n>Seq2\nAAA", ""), 'returncode': 0}
        process_mock.configure_mock(**attrs)
        mock_popen.return_value = process_mock

        mock_alignio_read.return_value = [
            SeqRecord(Seq("AAA"), id="Seq1"),
            SeqRecord(Seq("AAA"), id="Seq2")
        ]

        session_mock = MagicMock()
        mapping = UniProtPDBMapping(session_mock)

        # Ejecutar la función
        result = mapping.alinear_secuencias_mafft('AAA', 'AAA')

        # Verificar resultados y llamadas a los mocks
        self.assertEqual(result, ("AAA", "AAA", 100.0))
        mock_popen.assert_called()
        mock_alignio_read.assert_called()

    def test_calcular_porcentaje_identidad(self):
        session_mock = MagicMock()
        mapping = UniProtPDBMapping(session_mock)

        align = [
            SeqRecord(Seq("AAA"), id="Seq1"),
            SeqRecord(Seq("AAA"), id="Seq2")
        ]

        result = mapping.calcular_porcentaje_identidad(align)

        self.assertEqual(result, 100.0)

    def test_realizar_consulta_cadenas_iguales(self):
        # Configurar mock para el objeto session
        mock_session = MagicMock()
        mock_query = mock_session.query.return_value
        mock_query.join.return_value = mock_query
        mock_query.filter.return_value = mock_query
        mock_query.all.return_value = ['resultado simulado']

        mapping = UniProtPDBMapping(mock_session)

        # Ejecutar la función
        result = mapping.realizar_consulta_cadenas_iguales()

        # Verificar resultados y llamadas a los mocks
        self.assertEqual(result, ['resultado simulado'])
        mock_query.all.assert_called_once()
        mock_session.close.assert_called_once()

    @patch('concurrent.futures.ThreadPoolExecutor')
    @patch('protein_data_handler.alignment.UniProtPDBMapping.procesar_par')
    def test_volcar_datos_alineamiento(self, mock_procesar_par, mock_executor):
        mock_session = MagicMock()
        mapping = UniProtPDBMapping(mock_session)

        pares_simulados = [('par1', 'par2')]

        # Ejecutar la función
        mapping.volcar_datos_alineamiento(pares_simulados)

        # Verificar llamadas a los mocks
        mock_procesar_par.assert_called_with(('par1', 'par2'))
        mock_session.commit.assert_called_once()

    def test_procesar_par(self):
        mock_session = MagicMock()
        mapping = UniProtPDBMapping(mock_session)

        par_simulado = ('seq1', 'seq2', 'chain', 'pdb_id')
        mock_reference_id = 123
        mock_session.query.return_value.filter_by.return_value.one.return_value.id = mock_reference_id

        # Ejecutar la función
        mapping.procesar_par(par_simulado)

        # Verificar llamadas a los mocks
        mock_session.query.assert_called_with(PDBReference.id)
        mock_session.add.assert_called()


if __name__ == '__main__':
    unittest.main()
