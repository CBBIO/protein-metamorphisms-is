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
    """
    align_task_logger = logging.getLogger("align_task")
    try:
        print("Aligning structures using FATCAT...")

        pdb_chains_path = os.path.join(conf['data_directory'], 'models')

        # Correct dictionary-style access for 'alignment_entry'
        representative_name = f"{alignment_entry['subcluster_1_file_path']}"
        representative_structure_path = os.path.join(pdb_chains_path, representative_name)

        target_name = f"{alignment_entry['subcluster_2_file_path']}"
        target_structure_path = os.path.join(pdb_chains_path, target_name)

        temp_dir = tempfile.mkdtemp()

        # Convert the CIF files to PDB and save in temp directory
        representative_name_pdb = f"{representative_name}.pdb"
        representative_pdb_path = os.path.join(temp_dir, representative_name_pdb)
        cif_to_pdb(representative_structure_path, representative_pdb_path)

        target_name_pdb = f"{target_name}.pdb"
        target_pdb_path = os.path.join(temp_dir, target_name_pdb)
        cif_to_pdb(target_structure_path, target_pdb_path)

        # Check if files were created successfully and log their details
        if not os.path.exists(representative_pdb_path):
            align_task_logger.error(f"Representative PDB file not found: {representative_pdb_path}")
            return alignment_entry['subcluster_entry_1_id'], {'error_message': 'Representative PDB file not found'}

        if not os.path.exists(target_pdb_path):
            align_task_logger.error(f"Target PDB file not found: {target_pdb_path}")
            return alignment_entry['subcluster_entry_1_id'], {'error_message': 'Target PDB file not found'}

        # Build the FATCAT command
        # Build the FATCAT command
        command = [
            os.path.join(conf['binaries_path'], "FATCAT"),
            "-i", temp_dir,
            "-p1", representative_name_pdb,
            "-p2", target_name_pdb,
            "-b", "-q"
        ]

        # Log the full command for debugging
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Wait for the process to finish and capture the output and errors
        stdout, stderr = process.communicate()

        # Check if the command execution was successful
        if process.returncode != 1:
            print(f"Failed command: {' '.join(command)}")
            align_task_logger.error(f"FATCAT failed with return code {process.returncode} and stderr: {stderr}")
            return alignment_entry['subcluster_entry_1_id'], {
                'cluster_entry_id': alignment_entry['cluster_id'],
                'error_message': f"FATCAT failed with return code {process.returncode} and stderr: {stderr}"
            }

        # Parsing the output if the command succeeded
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

                'cluster_id': alignment_entry['cluster_id'],
                'subcluster_entry_1_id': alignment_entry['subcluster_entry_1_id'],
                'subcluster_entry_2_id': alignment_entry['subcluster_entry_2_id'],
                'fc_rms': rms,
                'fc_identity': identity,
                'fc_similarity': similarity,
                'fc_score': score,
                'fc_align_len': align_len
            }
        else:
            result = {
                'cluster_id': alignment_entry['cluster_id'],
                'subcluster_entry_1_id': alignment_entry['subcluster_entry_1_id'],
                'subcluster_entry_2_id': alignment_entry['subcluster_entry_2_id'],
                'error_message': stderr
            }


        return result

    except:
        pass