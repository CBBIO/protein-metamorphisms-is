from unittest import TestCase, mock
from unittest.mock import patch, MagicMock, mock_open

import pandas

from protein_metamorphisms_is.operations.cdhit import CDHit
from protein_metamorphisms_is.sql.model import PDBChains, Cluster


class TestCDHit(TestCase):
    def setUp(self):
        # Configuration for CDHit instance
        self.conf = {
            "DB_USERNAME": "usuario",
            "DB_PASSWORD": "clave",
            "DB_HOST": "localhost",
            "DB_PORT": 5432,
            "DB_NAME": "BioData",
            'fasta_path': './test.fasta',
            'cdhit_out_path': './test',
            'sequence_identity_threshold': 0.9,
            'alignment_coverage': 0.8,
            'memory_usage': 4000,
            'max_workers': 2,
            'most_representative_search': 1,
            'constants': 'path/to/constants.yaml'
        }

        with patch('protein_metamorphisms_is.information_system.base.extractor.setup_logger') as mock_logger, \
                patch('protein_metamorphisms_is.sql.base.database_manager.create_engine'), \
                patch('protein_metamorphisms_is.sql.base.database_manager.sessionmaker', return_value=MagicMock()), \
                patch('protein_metamorphisms_is.operations.base.operator.handle_structural_alignment_types'), \
                patch('protein_metamorphisms_is.operations.base.operator.handle_structural_complexity_levels'), \
                patch('yaml.safe_load'), \
                patch('builtins.open', mock.mock_open(read_data="constants: value")):
            self.cdhit = CDHit(self.conf)
            self.cdhit.logger = mock_logger

    @patch('protein_metamorphisms_is.operations.cdhit.CDHit.load_chains')
    @patch('protein_metamorphisms_is.operations.cdhit.CDHit.create_fasta')
    @patch('protein_metamorphisms_is.operations.cdhit.CDHit.cluster')
    def test_start_success(self, mock_cluster, mock_create_fasta, mock_load_chains):
        # Setup mock return values
        mock_load_chains.return_value = [MagicMock(), MagicMock()]

        # Execute the start method
        self.cdhit.start()

        # Verify that internal methods were called as expected
        mock_load_chains.assert_called_once()
        mock_create_fasta.assert_called_once_with(mock_load_chains.return_value)
        mock_cluster.assert_called_once()

    @patch('protein_metamorphisms_is.operations.cdhit.CDHit.load_chains')
    @patch('protein_metamorphisms_is.operations.cdhit.CDHit.create_fasta')
    @patch('protein_metamorphisms_is.operations.cdhit.CDHit.cluster')
    def test_start_exception(self, mock_cluster, mock_create_fasta, mock_load_chains):
        mock_load_chains.side_effect = Exception("Test exception")

        with self.assertRaises(Exception) as context:
            self.cdhit.start()

        self.assertEqual(str(context.exception), "Test exception")

        self.cdhit.logger.error.assert_called_with("Error during clustering process: Test exception")

    def test_load_chains_single_model(self):
        mock_chains = [mock.MagicMock(spec=PDBChains), mock.MagicMock(spec=PDBChains)]
        self.cdhit.session.query.return_value.filter.return_value.all.return_value = mock_chains

        chains = self.cdhit.load_chains()

        self.assertEqual(chains, mock_chains)
        self.cdhit.session.query.assert_called_once_with(PDBChains)

    def test_load_chains_multiple_models(self):
        self.cdhit.conf['allow_multiple_chain_models'] = True

        mock_chains = [mock.MagicMock(spec=PDBChains), mock.MagicMock(spec=PDBChains)]
        self.cdhit.session.query.return_value.all.return_value = mock_chains

        chains = self.cdhit.load_chains()

        self.assertEqual(chains, mock_chains)
        self.cdhit.session.query.assert_called_once_with(PDBChains)
        self.cdhit.session.query.return_value.filter.assert_not_called()

        self.cdhit.conf['allow_multiple_chain_models'] = False

    @patch('builtins.open', new_callable=mock_open)
    def test_create_fasta(self, mock_file_open):
        # Configurar los datos de prueba
        mock_chains = [PDBChains(id='chain1', sequence='AAA'), PDBChains(id='chain2', sequence='BBB')]

        # Ejecutar el método bajo prueba
        self.cdhit.create_fasta(mock_chains)

        # Verificar que se abrió el archivo con el path correcto y modo de escritura
        mock_file_open.assert_called_once_with(self.conf['fasta_path'], "w")

        # Preparar el contenido esperado del archivo FASTA
        expected_calls = [mock.call(f">{chain.id}\n{chain.sequence}\n") for chain in mock_chains]

        # Obtener el handle del archivo mockeado y verificar que se escribieron los datos esperados
        mock_file_handle = mock_file_open()
        mock_file_handle.write.assert_has_calls(expected_calls, any_order=True)

        # Verificar que se registró la acción en el logger
        self.cdhit.logger.info.assert_called_with(f"Writing protein chains to FASTA file at {self.conf['fasta_path']}")

    @patch('protein_metamorphisms_is.operations.cdhit.cd_hit')
    @patch('protein_metamorphisms_is.operations.cdhit.read_clstr')
    def test_cluster(self, mock_read_clstr, mock_cd_hit):
        # Configurar el mock para simular la lectura de los resultados del clúster
        mock_read_clstr.return_value = pandas.DataFrame([
            {"identifier": "chain1", "cluster": 1, "is_representative": True, "size": 100, "identity": 0.95},
            {"identifier": "chain2", "cluster": 1, "is_representative": False, "size": 90, "identity": 0.90}
        ])

        # Ejecutar el método bajo prueba
        self.cdhit.cluster()

        # Verificar que se llamó a cd_hit con los parámetros correctos
        mock_cd_hit.assert_called_once_with(
            i=self.conf['fasta_path'],
            o=self.conf['cdhit_out_path'],
            c=self.conf['sequence_identity_threshold'],
            d=0,
            sc=1,
            aL=self.conf['alignment_coverage'],
            M=self.conf['memory_usage'],
            T=self.conf['max_workers'],
            g=self.conf['most_representative_search']
        )

        # Verificar que se leyó el archivo de salida de CD-HIT
        mock_read_clstr.assert_called_once_with(f"{self.conf['cdhit_out_path']}.clstr")

        # Verificar que los resultados del clúster se almacenaron en la base de datos
        expected_calls = [
            mock.call.add(mock.ANY),  # mock.ANY se usa aquí como un comodín
            mock.call.commit()
        ]
        self.cdhit.session.assert_has_calls(expected_calls, any_order=True)

        # Verificar que se registraron las acciones en el logger
        expected_log_calls = [
            mock.call.info(f"Running CD-HIT on {self.conf['fasta_path']}"),
            mock.call.info(f"Reading CD-HIT output from {self.conf['cdhit_out_path']}.clstr"),
            mock.call.info("CD-HIT clustering data stored in the database")
        ]
        self.cdhit.logger.assert_has_calls(expected_log_calls)
