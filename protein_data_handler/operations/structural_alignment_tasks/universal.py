import logging
import os
import re
import subprocess


def align_task(alignment_entry, conf):
    """
    Performs the alignment task for a given pair of protein structures using TM-align.

    This function aligns a target protein structure with a representative structure using TM-align. It logs the
    process and handles any exceptions that occur during the alignment. The function is designed to be run asynchronously
    in a multiprocessing environment.

    Args:
        alignment_entry (dict): A dictionary containing the representative and target protein structures' information and paths.

    Returns:
        dict: A dictionary with the alignment result or error message, if any.
    """
    align_task_logger = logging.getLogger("align_task")
    try:
        align_task_logger.info("Aligning structures using TM-align...")
        pdb_chains_path = conf['pdb_chains_path']
        representative_name = f"{alignment_entry.rep_pdb_id}_{alignment_entry.rep_chains}_{alignment_entry.rep_model}"
        representative_structure_path = os.path.join(pdb_chains_path, f"{representative_name}.cif")
        target_name = f"{alignment_entry.pdb_id}_{alignment_entry.chains}_{alignment_entry.model}"
        target_structure_path = os.path.join(pdb_chains_path, f"{target_name}.cif")

        # Construye el comando para ejecutar TMalign
        command = ["TMalign", representative_structure_path, target_structure_path]
        print(command)
        # Ejecuta el comando
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Lee la salida y el error (si existe)
        stdout, stderr = process.communicate()

        if process.returncode == 0:
            # Procesa la salida si el comando fue exitoso
            # Expresiones regulares para buscar los valores deseados
            rmsd_pattern = re.compile(r"RMSD=\s*(\d+\.\d+)")
            seq_id_pattern = re.compile(r"Seq_ID=n_identical/n_aligned=\s*(\d+\.\d+)")
            tm_score_chain_1_pattern = re.compile(r"TM-score=\s*(\d+\.\d+)\s*\(if normalized by length of Chain_1\)")
            tm_score_chain_2_pattern = re.compile(r"TM-score=\s*(\d+\.\d+)\s*\(if normalized by length of Chain_2\)")

            rmsd = float(rmsd_pattern.search(stdout).group(1)) if rmsd_pattern.search(stdout) else None
            seq_id = float(seq_id_pattern.search(stdout).group(1)) if seq_id_pattern.search(stdout) else None
            tm_score_chain_1 = float(
                tm_score_chain_1_pattern.search(stdout).group(1)) if tm_score_chain_1_pattern.search(stdout) else None
            tm_score_chain_2 = float(
                tm_score_chain_2_pattern.search(stdout).group(1)) if tm_score_chain_2_pattern.search(stdout) else None

            result = {
                'cluster_entry_id': alignment_entry.cluster_id,
                'rmsd': rmsd,
                'seq_id': seq_id,
                'tm_score_chain_1': tm_score_chain_1,
                'tm_score_chain_2': tm_score_chain_2,
                'error_message': None
            }
            print(result)

        else:
            result = {
                'cluster_entry_id': alignment_entry.cluster_id,
                'error_message': stderr
            }

        align_task_logger.info("Alignment completed successfully.")
        return result
    except Exception as e:
        align_task_logger.error(f"Error during alignment task: {str(e)}")
        return {
            'cluster_entry_id': alignment_entry.cluster_id,
            'error_message': str(e)
        }
