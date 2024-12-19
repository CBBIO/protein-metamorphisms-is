import os
import tempfile
import subprocess
import re
import traceback

from protein_metamorphisms_is.helpers.parser.parser import cif_to_pdb


def align_task(alignment_entry, conf, logger):
    try:
        # Extract file paths and names from alignment_entry
        subcluster_1_path = alignment_entry['subcluster_1_file_path']
        subcluster_2_path = alignment_entry['subcluster_2_file_path']
        subcluster_1_name = os.path.basename(subcluster_1_path).replace('.cif', '')
        subcluster_2_name = os.path.basename(subcluster_2_path).replace('.cif', '')

        # Create a temporary directory for PDB files
        temp_dir = tempfile.mkdtemp()

        # Convert CIF to PDB
        subcluster_1_pdb_path = os.path.join(temp_dir, f"{subcluster_1_name}.pdb")
        cif_to_pdb(subcluster_1_path, subcluster_1_pdb_path)

        subcluster_2_pdb_path = os.path.join(temp_dir, f"{subcluster_2_name}.pdb")
        cif_to_pdb(subcluster_2_path, subcluster_2_pdb_path)

        # Check that the PDB files were created successfully
        if not os.path.exists(subcluster_1_pdb_path):
            raise FileNotFoundError(f"Failed to convert CIF to PDB for: {subcluster_1_path}")
        if not os.path.exists(subcluster_2_pdb_path):
            raise FileNotFoundError(f"Failed to convert CIF to PDB for: {subcluster_2_path}")

        # Prepare the FATCAT alignment command
        fatcat_command = [
            os.path.join(conf['binaries_path'], "FATCAT"),
            "-i", temp_dir,
            "-p1", f"{subcluster_1_name}.pdb",
            "-p2", f"{subcluster_2_name}.pdb",
            "-b", "-q"
        ]

        # Log and execute the FATCAT command
        logger.debug(f"Running FATCAT command: {' '.join(fatcat_command)}")
        process = subprocess.Popen(fatcat_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()

        # Log the outputs from the command
        logger.debug(f"FATCAT stdout: {stdout}")
        logger.debug(f"FATCAT stderr: {stderr}")

        # Check for execution success
        if process.returncode != 1:
            raise RuntimeError(f"FATCAT failed with return code {process.returncode}: {stderr}")

        # Parse alignment results
        rmsd_pattern = re.compile(r"opt-rmsd\s+(\d+\.\d+)")
        identity_pattern = re.compile(r"Identity\s+(\d+\.\d+)%")
        similarity_pattern = re.compile(r"Similarity\s+(\d+\.\d+)%")
        score_pattern = re.compile(r"Score\s+(\d+\.\d+)")
        align_len_pattern = re.compile(r"align-len\s+(\d+)")

        try:
            rmsd = float(rmsd_pattern.search(stdout).group(1)) if rmsd_pattern.search(stdout) else None
            identity = float(identity_pattern.search(stdout).group(1)) if identity_pattern.search(stdout) else None
            similarity = float(similarity_pattern.search(stdout).group(1)) if similarity_pattern.search(stdout) else None
            score = float(score_pattern.search(stdout).group(1)) if score_pattern.search(stdout) else None
            align_len = int(align_len_pattern.search(stdout).group(1)) if align_len_pattern.search(stdout) else None
        except Exception as parse_error:
            logger.error(f"Failed to parse FATCAT output: {parse_error}")
            logger.error(f"Raw output: {stdout}")
            raise

        # Return alignment result with cluster and subcluster IDs
        result = {
            'cluster_id': alignment_entry.get('cluster_id'),
            'subcluster_entry_1_id': alignment_entry.get('subcluster_entry_1_id'),
            'subcluster_entry_2_id': alignment_entry.get('subcluster_entry_2_id'),
            'fc_rms': rmsd,
            'fc_identity': identity,
            'fc_similarity': similarity,
            'fc_score': score,
            'fc_align_len': align_len
        }

        logger.info(f"Alignment result: {result}")
        return result

    except Exception as e:
        # Log the detailed error message with traceback
        error_message = f"Error during FATCAT alignment for cluster {alignment_entry.get('cluster_id')}: {str(e)}"
        logger.error(error_message)
        logger.error("Traceback:\n" + traceback.format_exc())
        return {
            'cluster_id': alignment_entry.get('cluster_id'),
            'subcluster_entry_1_id': alignment_entry.get('subcluster_entry_1_id'),
            'subcluster_entry_2_id': alignment_entry.get('subcluster_entry_2_id'),
            'error_message': error_message
        }
