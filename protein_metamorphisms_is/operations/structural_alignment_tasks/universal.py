import logging
import os
import re
import subprocess


def align_task(alignment_entry, conf):
    """
    Executes a protein structure alignment task using the US-align algorithm.

    This function aligns a target protein structure against a representative
    structure using US-align, an advanced algorithm for measuring structural
    similarity between protein structures. It runs an external US-align binary,
    captures its output, and extracts key metrics such as RMSD (Root Mean Square
    Deviation), sequence identity, and alignment scores.

    File paths for the representative and target structures are constructed based
    on their PDB IDs, chain identifiers, and model numbers. The US-align binary is
    executed with these structures as input, and its output is parsed using regular
    expressions to extract alignment metrics.

    If execution is successful, a dictionary containing these metrics is returned.
    If US-align encounters an error or if an exception occurs, the error is logged
    and an error object with the error message is returned.

    Args:
        alignment_entry (object): Contains data for the alignment task, including
            PDB IDs, chain identifiers, model numbers for both representative and
            target structures, and the cluster entry ID.
        conf (dict): Configuration settings, including paths to the directory where
            PDB chain files are stored and the directory containing the US-align
            binary.

    Returns:
        tuple: Contains the queue entry ID of the alignment task and either a result
            dictionary with keys for 'cluster_entry_id', 'us_rms', 'us_seq_id',
            'us_score', or an error object with 'cluster_entry_id' and
            'error_message'.

    Raises:
        Exception: Any exceptions during the process are captured, logged, and
            an error object with the error message is returned.

    Example:
        >>> alignment_entry = {
                'rep_pdb_id': '1A2B',
                'rep_chains': 'A',
                'rep_model': 0,
                'pdb_id': '2B3C',
                'chains': 'B',
                'model': 0,
                'cluster_id': 123,
                'queue_entry_id': 456
            }
        >>> conf = {
                'pdb_chains_path': '/path/to/pdb_chains',
                'binaries_path': '/path/to/binaries'
            }
        >>> align_task(alignment_entry, conf)
        (456, {'cluster_entry_id': 123, 'us_rms': 0.5, 'us_seq_id': 0.8, 'us_score': 0.95})
    """
    align_task_logger = logging.getLogger("align_task")
    try:
        align_task_logger.info("Aligning structures using US-align...")
        pdb_chains_path = conf['pdb_chains_path']
        representative_name = f"{alignment_entry.rep_pdb_id}_{alignment_entry.rep_chains}_{alignment_entry.rep_model}"
        representative_structure_path = os.path.join(pdb_chains_path, f"{representative_name}.cif")
        target_name = f"{alignment_entry.pdb_id}_{alignment_entry.chains}_{alignment_entry.model}"
        target_structure_path = os.path.join(pdb_chains_path, f"{target_name}.cif")

        command = [os.path.join(conf['binaries_path'], "USalign"), representative_structure_path, target_structure_path]

        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Lee la salida y el error (si existe)
        stdout, stderr = process.communicate()
        if process.returncode == 0:
            # Procesa la salida si el comando fue exitoso
            # Expresiones regulares para buscar los valores deseados
            rmsd_pattern = re.compile(r"RMSD=\s*(\d+\.\d+)")
            seq_id_pattern = re.compile(r"Seq_ID=n_identical/n_aligned=\s*(\d+\.\d+)")
            tm_score_chain_1_pattern = re.compile(r"TM-score=\s*(\d+\.\d+)\s*\(normalized by length of Structure_1")
            tm_score_chain_2_pattern = re.compile(r"TM-score=\s*(\d+\.\d+)\s*\(normalized by length of Structure_2")

            rmsd = float(rmsd_pattern.search(stdout).group(1)) if rmsd_pattern.search(stdout) else None
            seq_id = float(seq_id_pattern.search(stdout).group(1)) if seq_id_pattern.search(stdout) else None
            tm_score_chain_1 = float(
                tm_score_chain_1_pattern.search(stdout).group(1)) if tm_score_chain_1_pattern.search(stdout) else None
            tm_score_chain_2 = float(
                tm_score_chain_2_pattern.search(stdout).group(1)) if tm_score_chain_2_pattern.search(stdout) else None

            result = {
                'cluster_entry_id': alignment_entry.cluster_id,
                'tm_rms': rmsd,
                'tm_seq_id': seq_id,
                'tm_score_chain_1': tm_score_chain_1,
                'tm_score_chain_2': tm_score_chain_2
            }

        else:
            result = {
                'cluster_entry_id': alignment_entry.cluster_id,
                'error_message': stderr
            }

        align_task_logger.info("Alignment completed successfully.")
        return alignment_entry.queue_entry_id, result
    except Exception as e:
        align_task_logger.error(f"Error during alignment task: {str(e)}")
        error_object = {'error_message': e}
        return alignment_entry.queue_entry_id, error_object
