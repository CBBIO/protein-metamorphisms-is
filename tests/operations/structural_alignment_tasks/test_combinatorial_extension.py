from unittest.mock import patch, MagicMock

from protein_metamorphisms_is.operation.structural_alignment_tasks.combinatorial_extension import align_task


class AlignmentEntry:
    def __init__(self, cluster_id, rep_pdb_id, rep_chains, rep_model, pdb_id, chains, model,
                 queue_entry_id):
        self.cluster_id = cluster_id
        self.rep_pdb_id = rep_pdb_id
        self.rep_chains = rep_chains
        self.rep_model = rep_model
        self.pdb_id = pdb_id
        self.chains = chains
        self.model = model
        self.queue_entry_id = queue_entry_id


def test_align_task_successful():
    with patch(
            'protein_metamorphisms_is.operations.structural_alignment_tasks.combinatorial_extension.os.path.join') as mock_join:
        with patch(
                'protein_metamorphisms_is.operations.structural_alignment_tasks.combinatorial_extension.MMCIFParser') as mock_parser:
            # Configura mocks para simular el comportamiento esperado
            mock_join.return_value = '/fake/path/to/file.cif'
            mock_structure = MagicMock()
            mock_parser.return_value.get_structure.return_value = mock_structure

            # Suponiendo que CEAligner y su método `align` pueden ser mockeados adecuadamente
            with patch(
                    'protein_metamorphisms_is.operations.structural_alignment_tasks.combinatorial_extension.CEAligner') as mock_aligner:
                mock_aligner.return_value.rms = 0.5

                # Creando una instancia de la clase con los mismos valores que el diccionario original
                alignment_entry = AlignmentEntry(
                    cluster_id=123,
                    rep_pdb_id='1A2B',
                    rep_chains='A',
                    rep_model=0,
                    pdb_id='2B3C',
                    chains='B',
                    model=0,
                    queue_entry_id=456
                )

                conf = {'pdb_chains_path': '/fake/path'}

                # Llamada a la función bajo prueba
                queue_id, result = align_task(alignment_entry, conf)

                # Verificaciones
                assert queue_id == 456
                assert result == {'cluster_entry_id': 123, 'ce_rms': 0.5}


def test_align_task_exception():
    with patch(
            'protein_metamorphisms_is.operations.structural_alignment_tasks.combinatorial_extension.os.path.join') as mock_join, \
            patch(
                'protein_metamorphisms_is.operations.structural_alignment_tasks.combinatorial_extension.MMCIFParser') as mock_parser, \
            patch(
                'protein_metamorphisms_is.operations.structural_alignment_tasks.combinatorial_extension.CEAligner') as mock_aligner:
        mock_join.return_value = '/fake/path/to/file.cif'
        mock_aligner.side_effect = Exception("Test error")
        # Datos de entrada simulados usando la clase AlignmentEntry
        alignment_entry = AlignmentEntry(
            cluster_id=123,
            rep_pdb_id='1A2B',
            rep_chains='A',
            rep_model=0,
            pdb_id='2B3C',
            chains='B',
            model=0,
            queue_entry_id=456
        )
        conf = {'pdb_chains_path': '/fake/path'}
        # Llamada a la función bajo prueba
        queue_id, result = align_task(alignment_entry, conf)
        # Verificaciones
        assert queue_id == 456
        assert 'error_message' in result
        assert result['cluster_entry_id'] == 123
        assert isinstance(result['error_message'], Exception)
