import os
import re
import subprocess
import traceback


def align_task(alignment_entry, conf, logger):
    """
    Executes a protein structure alignment task using the US-align algorithm.
    """
    try:
        logger.info("Aligning structures using US-align...")

        pdb_chains_path = os.path.join(conf['data_directory'], 'models')

        # Access the subcluster file paths
        representative_name = alignment_entry.get('subcluster_1_file_path')
        target_name = alignment_entry.get('subcluster_2_file_path')

        if not representative_name or not target_name:
            raise ValueError("Missing file paths for representative or target structures.")

        representative_structure_path = os.path.join(pdb_chains_path, representative_name)
        target_structure_path = os.path.join(pdb_chains_path, target_name)

        # Construct the US-align command
        command = [os.path.join(conf['binaries_path'], "USalign"), representative_structure_path, target_structure_path]
        logger.info(f"Running US-align command: {' '.join(command)}")

        # Execute the US-align command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()

        # Check the result of the process
        if process.returncode == 0:
            # Regular expressions to extract the values from the US-align output
            rmsd_pattern = re.compile(r"RMSD=\s*(\d+\.\d+)")
            seq_id_pattern = re.compile(r"Seq_ID=n_identical/n_aligned=\s*(\d+\.\d+)")
            tm_score_chain_1_pattern = re.compile(r"TM-score=\s*(\d+\.\d+)\s*\(normalized by length of Structure_1")
            tm_score_chain_2_pattern = re.compile(r"TM-score=\s*(\d+\.\d+)\s*\(normalized by length of Structure_2")

            rmsd = float(rmsd_pattern.search(stdout).group(1)) if rmsd_pattern.search(stdout) else None
            seq_id = float(seq_id_pattern.search(stdout).group(1)) if seq_id_pattern.search(stdout) else None
            tm_score_chain_1 = float(tm_score_chain_1_pattern.search(stdout).group(1)) if tm_score_chain_1_pattern.search(stdout) else None
            tm_score_chain_2 = float(tm_score_chain_2_pattern.search(stdout).group(1)) if tm_score_chain_2_pattern.search(stdout) else None

            # Return the results as a dictionary
            result = {
                'cluster_id': alignment_entry.get('cluster_id'),
                'subcluster_entry_1_id': alignment_entry.get('subcluster_entry_1_id'),
                'subcluster_entry_2_id': alignment_entry.get('subcluster_entry_2_id'),
                'tm_rms': rmsd,
                'tm_seq_id': seq_id,
                'tm_score_chain_1': tm_score_chain_1,
                'tm_score_chain_2': tm_score_chain_2
            }

            logger.info("Alignment completed successfully.")
        else:
            # Handle process error
            error_message = f"US-align failed with return code {process.returncode} and stderr: {stderr.strip()}"
            logger.error(error_message)
            result = {
                'cluster_id': alignment_entry.get('cluster_id'),
                'subcluster_entry_1_id': alignment_entry.get('subcluster_entry_1_id'),
                'subcluster_entry_2_id': alignment_entry.get('subcluster_entry_2_id'),
                'error_message': error_message
            }

        return result

    except Exception as e:
        # Log detailed error message with traceback
        error_message = f"Error during US-alignment for cluster {alignment_entry.get('cluster_id')}: {str(e)}"
        logger.error(error_message)
        logger.error("Traceback:\n" + traceback.format_exc())
        return {
            'cluster_id': alignment_entry.get('cluster_id'),
            'subcluster_entry_1_id': alignment_entry.get('subcluster_entry_1_id'),
            'subcluster_entry_2_id': alignment_entry.get('subcluster_entry_2_id'),
            'error_message': error_message
        }
