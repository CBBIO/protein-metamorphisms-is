import logging
import os

from Bio.PDB import MMCIFParser, CEAligner


def align_task(alignment_entry, conf):
    """
    Performs the alignment task for a given pair of protein structures.

    This function aligns a target protein structure with a representative structure using the Combinatorial Extension (CE) algorithm. It is designed to be run as a part of a parallel processing system for structural alignments, managing the alignment of multiple protein structures efficiently.

    The function retrieves the paths for the representative and target protein structures based on their PDB IDs, chain identifiers, and model numbers from the configuration settings. It then parses the structures using Biopython's MMCIFParser, sets the representative structure as a reference, and performs the alignment against the target structure. The root-mean-square deviation (RMSD) of the alignment is calculated and returned along with the cluster entry ID as a measure of the alignment's quality.

    In case of any exceptions during the alignment process, an error message is logged, and an error object is returned, containing the cluster entry ID and the error message.

    Args:
        alignment_entry (object): An object containing the alignment task's data, including PDB IDs, chain identifiers, and model numbers for the representative and target structures, as well as the cluster entry ID.
        conf (dict): A dictionary containing configuration settings, specifically the path to the directory where PDB chain files are stored.

    Returns:
        tuple: A tuple containing the queue entry ID of the alignment task and either the alignment result (a dictionary with the cluster entry ID and the CE RMSD value) or an error object (a dictionary with the cluster entry ID and an error message).

    Raises:
        Exception: Captures and logs any exceptions that occur during the alignment process, returning an error object with details.

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
        >>> conf = {'pdb_chains_path': '/path/to/pdb_chains'}
        >>> align_task(alignment_entry, conf)
        (456, {'cluster_entry_id': 123, 'ce_rms': 0.5})
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
            'ce_rms': rms
        }

        align_task_logger.info("Alignment completed successfully.")
        return alignment_entry.queue_entry_id, result
    except Exception as e:
        error_object = {'cluster_entry_id': alignment_entry.cluster_id, 'error_message': e}
        return alignment_entry.queue_entry_id, error_object
