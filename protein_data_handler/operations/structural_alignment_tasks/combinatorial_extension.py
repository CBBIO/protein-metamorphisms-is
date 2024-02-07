import logging
import os

from Bio.PDB import MMCIFParser, CEAligner


def align_task(alignment_entry,conf):
    """
    Performs the alignment task for a given pair of protein structures.

    This function aligns a target protein structure with a representative structure using the CE algorithm. It logs the
    process and handles any exceptions that occur during the alignment. The function is designed to be run asynchronously
    in a multiprocessing environment.

    Args:
        alignment_entry (dict): A dictionary containing the representative and target protein structures' information and paths.

    Returns:
        dict: A dictionary with the alignment result or error message, if any.
    """
    align_task_logger = logging.getLogger("align_task")
    try:
        align_task_logger.info("Aligning structures...")
        pdb_chains_path = conf['pdb_chains_path']

        representative_name = f"{alignment_entry.rep_pdb_id}_{alignment_entry.rep_chains}_{alignment_entry.rep_model}"
        representative_structure_path = os.path.join(pdb_chains_path, f"{representative_name}.cif")
        target_name = f"{alignment_entry.pdb_id}_{alignment_entry.chains}_{alignment_entry.model}"
        target_structure_path = os.path.join(pdb_chains_path, f"{target_name}.cif")

        parser = MMCIFParser()
        representative_structure = parser.get_structure(representative_name, representative_structure_path)
        aligner = CEAligner()
        aligner.set_reference(representative_structure)
        target_structure = parser.get_structure(target_name, target_structure_path)

        aligner.align(target_structure)
        rms = aligner.rms

        result = {
            'cluster_entry_id': alignment_entry.cluster_id,
            'rms': rms
        }
        print(result)

        align_task_logger.info("Alignment completed successfully.")
        return result
    except Exception as e:
        error_object = {'cluster_entry_id': alignment_entry.cluster_id, 'error_message': e}
        return error_object