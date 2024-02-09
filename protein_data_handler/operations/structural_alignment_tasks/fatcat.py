import logging
import os
import re
import subprocess
import tempfile
import traceback

from protein_data_handler.helpers.parser.parser import cif_to_pdb


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

        temp_dir = tempfile.mkdtemp()

        # Convertir los archivos CIF a PDB y guardarlos en el directorio temporal
        representative_name_pdb = f"{representative_name}.pdb"
        representative_pdb_path = os.path.join(temp_dir, representative_name_pdb)
        cif_to_pdb(representative_structure_path, representative_pdb_path)
        target_name_pdb = f"{target_name}.pdb"
        target_pdb_path = os.path.join(temp_dir, target_name_pdb)
        cif_to_pdb(target_structure_path, target_pdb_path)

        # Construye el comando para ejecutar TMalign
        command = [
            os.path.join(conf['binaries_path'], "FATCAT"),
            "-i", temp_dir,
            "-p1", representative_name_pdb,
            "-p2", target_name_pdb,
            "-b", "-q"
        ]

        # Ejecuta el comando
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Lee la salida y el error (si existe)
        stdout, stderr = process.communicate()
        if process.returncode == 1:
            rmsd_pattern = re.compile(r"opt-rmsd\s+(\d+\.\d+)")
            identity_pattern = re.compile(r"Identity\s+(\d+\.\d+)%")
            similarity_pattern = re.compile(r"Similarity\s+(\d+\.\d+)%")
            score_pattern = re.compile(r"Score\s+(\d+\.\d+)")
            align_len_pattern = re.compile(r"align-len\s+(\d+)")

            rms = float(rmsd_pattern.search(stdout).group(1)) if rmsd_pattern.search(stdout) else None
            identity = float(identity_pattern.search(stdout).group(1)) if identity_pattern.search(stdout) else None
            similarity = float(similarity_pattern.search(stdout).group(1)) if similarity_pattern.search(
                stdout) else None
            score = float(score_pattern.search(stdout).group(1)) if score_pattern.search(stdout) else None
            align_len = int(align_len_pattern.search(stdout).group(1)) if align_len_pattern.search(stdout) else None

            result = {
                'cluster_entry_id': alignment_entry.cluster_id,
                'fc_rms': rms,
                'fc_identity': identity,
                'fc_similarity': similarity,
                'fc_score': score,
                'fc_align_len': align_len
            }

        else:
            result = {
                'cluster_entry_id': alignment_entry.cluster_id,
                'error_message': stderr
            }

        align_task_logger.info("Alignment completed successfully.")
        return alignment_entry.queue_entry_id, result
    except Exception as e:
        traceback_info = traceback.format_exc()
        align_task_logger.error(f"Error during alignment task: {str(e)} Traceback:\n{traceback_info}")

        error_object = {'error_message': e}
        return alignment_entry.queue_entry_id, error_object
