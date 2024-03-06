import pytest
from unittest.mock import patch, MagicMock
from protein_metamorphisms_is.operations.structural_alignment_tasks.universal import align_task

@patch('subprocess.Popen')
def test_align_task_success(mock_popen):
    # Configura la salida simulada de USalign
    mock_process = MagicMock()
    stdout_data = """
    RMSD= 0.5
    Seq_ID=n_identical/n_aligned= 0.8
    TM-score= 0.95 (normalized by length of Structure_1)
    TM-score= 0.95 (normalized by length of Structure_2)
    """
    mock_process.communicate.return_value = (stdout_data, '')
    mock_process.returncode = 0
    mock_popen.return_value = mock_process

    # Configuración simulada de entrada
    alignment_entry = MagicMock()
    alignment_entry.cluster_id = 123
    alignment_entry.queue_entry_id = 456
    conf = {'pdb_chains_path': '/fake/path', 'binaries_path': '/fake/binaries'}

    queue_id, result = align_task(alignment_entry, conf)

    assert queue_id == 456
    assert result['cluster_entry_id'] == 123
    assert result['tm_rms'] == 0.5
    assert result['tm_seq_id'] == 0.8
    assert 'tm_score_chain_1' in result
    assert 'tm_score_chain_2' in result


@patch('subprocess.Popen')
def test_align_task_usalign_error(mock_popen):
    mock_process = MagicMock()
    mock_process.communicate.return_value = ('', 'Error message from USalign')
    mock_process.returncode = 1  # Un código de error
    mock_popen.return_value = mock_process

    alignment_entry = MagicMock()
    alignment_entry.cluster_id = 123
    alignment_entry.queue_entry_id = 456
    conf = {'pdb_chains_path': '/fake/path', 'binaries_path': '/fake/binaries'}

    queue_id, result = align_task(alignment_entry, conf)

    assert queue_id == 456
    assert 'error_message' in result


@patch('subprocess.Popen')
def test_align_task_exception(mock_popen):
    mock_popen.side_effect = Exception("Simulated exception")

    alignment_entry = MagicMock()
    alignment_entry.cluster_id = 123
    alignment_entry.queue_entry_id = 456
    conf = {'pdb_chains_path': '/fake/path', 'binaries_path': '/fake/binaries'}

    queue_id, result = align_task(alignment_entry, conf)

    assert queue_id == 456
    assert 'error_message' in result
