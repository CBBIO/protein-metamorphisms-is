import os
import traceback
from Bio.PDB import MMCIFParser, CEAligner
from Bio.PDB.PDBExceptions import PDBException


def align_task(alignment_entry, conf, logger):
    """
    Performs the alignment task for a given pair of protein structures using the
    Combinatorial Extension (CE) algorithm.
    """
    try:
        logger.info("Starting alignment using Combinatorial Extension (CE)...")

        # Build the paths for the PDB files
        pdb_chains_path = os.path.join(conf['data_directory'], 'models')

        representative_name = alignment_entry.get('subcluster_1_file_path')
        target_name = alignment_entry.get('subcluster_2_file_path')

        if not representative_name or not target_name:
            raise ValueError("Missing file paths for representative or target structures.")

        representative_structure_path = os.path.join(pdb_chains_path, representative_name)
        target_structure_path = os.path.join(pdb_chains_path, target_name)

        logger.info(f"Representative structure path: {representative_structure_path}")
        logger.info(f"Target structure path: {target_structure_path}")

        # Parse the structures
        parser = MMCIFParser()
        representative_structure = parser.get_structure(representative_name, representative_structure_path)
        target_structure = parser.get_structure(target_name, target_structure_path)

        # Perform the alignment using CEAligner
        aligner = CEAligner()

        try:
            aligner.set_reference(representative_structure)
            aligner.align(target_structure)

        except PDBException as e:
            # Check if the error is due to too few atoms in the reference structure
            if "Too few atoms in the reference structure" in str(e):
                logger.warning(f"{str(e)}. Reducing the window_size parameter.")
                # Here you would typically reduce the window_size; for example:
                aligner.window_size = 5  # Ajusta este valor según sea necesario
                # Retry alignment after adjusting the window size
                aligner.set_reference(representative_structure)
                aligner.align(target_structure)
            else:
                raise  # Re-raise if the exception is not related to too few atoms

        # Extract the RMS value from the alignment
        rms = aligner.rms

        # Build the result dictionary, propagating the cluster and subcluster IDs
        result = {
            'cluster_id': alignment_entry.get('cluster_id'),
            'subcluster_entry_1_id': alignment_entry.get('subcluster_entry_1_id'),
            'subcluster_entry_2_id': alignment_entry.get('subcluster_entry_2_id'),
            'ce_rms': rms
        }

        logger.info("Alignment completed successfully.")
        return result

    except Exception as e:
        # Log the error with traceback
        error_message = f"Error during CE-alignment for cluster {alignment_entry.get('cluster_id')}: {str(e)}"
        logger.error(error_message)
        logger.error("Traceback:\n" + traceback.format_exc())
        return {
            'cluster_id': alignment_entry.get('cluster_id'),
            'subcluster_entry_1_id': alignment_entry.get('subcluster_entry_1_id'),
            'subcluster_entry_2_id': alignment_entry.get('subcluster_entry_2_id'),
            'error_message': error_message
        }
