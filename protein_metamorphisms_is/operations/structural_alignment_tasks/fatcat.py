import logging
import os
import re
import subprocess
import tempfile
import traceback

from protein_metamorphisms_is.helpers.parser.parser import cif_to_pdb


def align_task(alignment_entry, conf):
    """
        Executes a protein structure alignment task using FATCAT, preceded by CIF to PDB conversion.

        This function aligns a target protein structure against a representative structure using the FATCAT algorithm,
        renowned for its flexibility in handling protein comparisons. Prior to alignment, the function converts CIF files
        to PDB format using a helper function, ensuring compatibility with FATCAT. The alignment process extracts several
        key metrics from the output, including RMSD, sequence identity, similarity, score, and alignment length.

        The function constructs paths for the representative and target structures, converts them to PDB format, and
        executes the FATCAT binary with these structures as input. It captures and parses the output to extract alignment
        metrics.

        In case of successful execution, it returns these metrics encapsulated in a result dictionary. If FATCAT encounters
        an error or if an exception occurs during the process, the function logs detailed error information including a
        traceback and returns an error object containing the error message.

        Args:
            alignment_entry (object): An object containing data for the alignment task, including PDB IDs, chain identifiers,
                                       and model numbers for both the representative and target structures, as well as the
                                       cluster entry ID.
            conf (dict): A dictionary containing configuration settings, specifically the paths to the directories where
                         PDB chain files are stored and where the FATCAT binary is located.

        Returns:
            tuple: A tuple containing the queue entry ID of the alignment task and either a result dictionary (with keys for
                   'cluster_entry_id', 'fc_rms', 'fc_identity', 'fc_similarity', 'fc_score', 'fc_align_len') or an error
                   object (with 'cluster_entry_id' and 'error_message').

        Raises:
            Exception: Logs any exceptions that occur during the process, including a traceback, returning an error object
                       with the error message.

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
            (456, {'cluster_entry_id': 123, 'fc_rms': 0.5, 'fc_identity': 99.0, 'fc_similarity': 98.5,
                   'fc_score': 150.0, 'fc_align_len': 250})
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
