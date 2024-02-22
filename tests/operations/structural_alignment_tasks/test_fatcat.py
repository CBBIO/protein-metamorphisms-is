import pytest
from unittest.mock import patch, MagicMock
from protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat import align_task

# Prueba de caso de éxito
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.os.path.join')
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.tempfile.mkdtemp')
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.cif_to_pdb')
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.subprocess.Popen')
def test_align_task_successful(mock_popen, mock_cif_to_pdb, mock_mkdtemp, mock_join):
    mock_mkdtemp.return_value = '/fake/temp_dir'
    mock_join.side_effect = lambda *args: '/fake/' + '/'.join(args)

    # Configura el mock de Popen para simular la salida de FATCAT
    process_mock = MagicMock()
    process_mock.communicate.return_value = ('output', '')
    process_mock.returncode = 1
    mock_popen.return_value = process_mock

    alignment_entry = MagicMock()
    alignment_entry.cluster_id = 123
    alignment_entry.queue_entry_id = 456
    alignment_entry.rep_pdb_id = '1A2B'
    alignment_entry.rep_chains = 'A'
    alignment_entry.rep_model = 0
    alignment_entry.pdb_id = '2B3C'
    alignment_entry.chains = 'B'
    alignment_entry.model = 0

    conf = {
        'pdb_chains_path': '/path/to/pdb_chains',
        'binaries_path': '/path/to/binaries'
    }

    queue_id, result = align_task(alignment_entry, conf)

    assert queue_id == 456
    assert 'fc_rms' in result
    assert 'fc_identity' in result
    assert 'fc_similarity' in result
    assert 'fc_score' in result
    assert 'fc_align_len' in result

# Prueba de caso de error
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.os.path.join')
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.tempfile.mkdtemp')
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.cif_to_pdb')
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.subprocess.Popen')
def test_align_task_error(mock_popen, mock_cif_to_pdb, mock_mkdtemp, mock_join):
    mock_mkdtemp.return_value = '/fake/temp_dir'
    mock_join.side_effect = lambda *args: '/fake/' + '/'.join(args)

    # Configura el mock de Popen para simular un error de FATCAT
    process_mock = MagicMock()
    process_mock.communicate.return_value = ('', 'error message')
    process_mock.returncode = 0  # Suponiendo que 0 indica un error en este contexto
    mock_popen.return_value = process_mock

    alignment_entry = MagicMock()
    alignment_entry.cluster_id = 123
    alignment_entry.queue_entry_id = 456

    conf = {
        'pdb_chains_path': '/path/to/pdb_chains',
        'binaries_path': '/path/to/binaries'
    }

    queue_id, result = align_task(alignment_entry, conf)

    assert queue_id == 456
    assert 'error_message' in result
    assert result['error_message'] == 'error message'


# Prueba que simula una excepción durante la tarea de alineación
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.os.path.join')
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.tempfile.mkdtemp')
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.cif_to_pdb')
@patch('protein_metamorphisms_is.operations.structural_alignment_tasks.fatcat.subprocess.Popen')
def test_align_task_exception(mock_popen, mock_cif_to_pdb, mock_mkdtemp, mock_join):
    mock_mkdtemp.return_value = '/fake/temp_dir'
    mock_join.side_effect = lambda *args: '/fake/' + '/'.join(args)

    # Configura mock_cif_to_pdb o mock_popen para lanzar una excepción
    mock_cif_to_pdb.side_effect = Exception("Simulated conversion error")

    alignment_entry = MagicMock()
    alignment_entry.cluster_id = 123
    alignment_entry.queue_entry_id = 456
    alignment_entry.rep_pdb_id = '1A2B'
    alignment_entry.rep_chains = 'A'
    alignment_entry.rep_model = 0
    alignment_entry.pdb_id = '2B3C'
    alignment_entry.chains = 'B'
    alignment_entry.model = 0

    conf = {
        'pdb_chains_path': '/path/to/pdb_chains',
        'binaries_path': '/path/to/binaries'
    }

    queue_id, result = align_task(alignment_entry, conf)

    # Verifica que la función maneja correctamente la excepción
    assert queue_id == 456
    assert 'error_message' in result
    # Dependiendo de cómo esté implementado tu manejo de errores, es posible que necesites ajustar esta aserción
    # para verificar el contenido específico del mensaje de error.
    assert isinstance(result['error_message'], Exception)