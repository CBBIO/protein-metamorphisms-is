import os
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio.PDB import PDBParser, CEAligner
from sqlalchemy.orm import sessionmaker

from protein_data_handler.operations.base.bioinfo_operator import BioinfoOperatorBase
from protein_data_handler.sql.model import PDBChains, Cluster, PDBReference, CEAlignResults


class CEAlign(BioinfoOperatorBase):
    """
    This class performs structural alignment of protein data using the Combinatorial Extension (CE) algorithm.
    CE algorithm is a popular method for protein structure alignment, known for its effectiveness in identifying
    resemblances between proteins that share very little sequence similarity. It works by creating an optimal
    alignment between two protein structures based on their backbone atom positions.
    """

    def __init__(self, conf):
        """
        Initializes the CEAlign object with configuration settings.

        Args:
            conf (dict): Configuration parameters, including database connections and operational settings.
        """
        super().__init__(conf, session_required=True)
        self.logger.info("CEAlign instance created with configuration.")

    def start(self):
        """
        Begins the structural alignment process. This method orchestrates the entire alignment workflow,
        starting from loading the cluster representatives, performing CE alignment, and handling any exceptions
        that occur during the process.
        """
        try:
            self.logger.info("Starting structural alignment process.")
            cluster_entries = self.load_clusters()

            alignment_map = map_representatives_to_targets(cluster_entries)

            self.logger.info(f"{len(alignment_map)} clusters avaliable")

            self.ce_align(alignment_map)

        except Exception as e:
            self.logger.error(f"Error during structural alignment process: {e}")
            raise

    def load_clusters(self):
        """
        Loads cluster representatives from the database. These representatives are typically selected protein
        structures that serve as a reference for their respective clusters in protein structural analysis.

        Returns:
            list: A list of cluster representative objects.
        """
        self.logger.info("Loading cluster entries from database.")
        clusters = self.session.query(
            Cluster.id,
            Cluster.cluster_id,
            Cluster.is_representative,
            PDBChains.chains,
            PDBReference.pdb_id
        ).join(
            PDBChains, Cluster.pdb_chain_id == PDBChains.id
        ).join(
            PDBReference, PDBChains.pdb_reference_id == PDBReference.id
        ).all()
        self.logger.info(f"Loaded {len(clusters)} cluster entries.")

        results = []
        for row in clusters:
            row_dict = {
                "id": row[0],
                "cluster_id": row[1],
                "is_representative": row[2],
                "chains": row[3],
                "pdb_id": row[4]
            }
            results.append(row_dict)
        return results

    def ce_align(self, alignment_map):

        num_workers = self.conf.get('max_workers', 4)
        self.logger.info(f"Performing CE alignment with {num_workers} workers.")
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(self.align_task, alignment_entry) for alignment_entry in alignment_map]

            # for future in as_completed(futures):
            #     try:
            #         alignment_results = future.result()
            #         self.logger.info("Alignment task completed successfully.")
            #     except Exception as e:
            #         self.logger.error(f"Error in alignment task: {e}", exc_info=True)

    def align_task(self, alignment_entry):
        """
        Task to align a single representative against all targets in its cluster. This function is where
        the actual structural alignment takes place using the CE algorithm. It includes loading the protein
        structures from PDB files, aligning them, and computing the RMS (Root Mean Square) deviation,
        which quantifies the similarity between the structures.

        Args:
            alignment_entry (tuple): Tuple containing the cluster_id and a dictionary with representative and targets.
        """
        representative = alignment_entry['representative']
        target = alignment_entry['target']
        cluster_id = representative['cluster_id']

        pdb_chains_path = self.conf.get('pdb_chains_path', './chains')
        parser = PDBParser()

        representative_name = f"{representative['pdb_id']}_{representative['chains']}"
        representative_structure_path = os.path.join(pdb_chains_path, f"{representative_name}.pdb")

        self.logger.info(f"Aligning representative structure: {representative_name} with cluster id {cluster_id}")

        try:
            representative_structure = parser.get_structure(representative_name, representative_structure_path)
            aligner = CEAligner()
            aligner.set_reference(representative_structure)


            target_name = f"{target['pdb_id']}_{target['chains']}"
            target_structure_path = os.path.join(pdb_chains_path, f"{target_name}.pdb")
            target_structure = parser.get_structure(target_name, target_structure_path)

            aligner.align(target_structure)
            self.logger.info(f'{cluster_id}')
            rms = aligner.rms

            result = CEAlignResults(cluster_entry_id=target['id'],
                                            rms=rms)  # Assuming 'id' is a field in target
            self.session.add(result)
            self.session.commit()
        except Exception as e:
            self.logger.error(
                f"Error in processing target {target['id']}: {e}")  # Assuming 'id' is a field in target


        except Exception as e:
            self.logger.error(
                f"Error in align_task for representative {representative['id']}: {e}")  # Assuming 'id' is a field in representative
            self.session.rollback()
        finally:
            self.session.close()

        self.logger.info(f"Alignment task for representative {representative_name} completed.")


def get_cluster_pdb_chain_details(self, cluster_entry_id):
    """
    Retrieves PDB chain details for a given cluster entry. This method queries the database to fetch
    information about the protein structure, such as its PDB ID and chain identifier. This information
    is crucial for locating the correct PDB file and understanding the context of the alignment.

    Args:
        cluster_entry_id (int): ID of the cluster entry.

    Returns:
        tuple: Details of the PDB chain associated with the cluster entry.
    """
    result = self.session.query(
        Cluster.id.label("cluster_id"),
        PDBChains.chains.label("chain"),
        PDBReference.pdb_id.label("pdb_id")
    ).join(
        PDBChains, Cluster.pdb_chain_id == PDBChains.id
    ).join(
        PDBReference, PDBChains.pdb_reference_id == PDBReference.id
    ).filter(
        Cluster.id == cluster_entry_id
    ).first()

    if result is None:
        self.logger.warning(f"No PDB chain details found for cluster entry {cluster_entry_id}.")
    else:
        self.logger.info(f"PDB chain details retrieved for cluster entry {cluster_entry_id}.")
    return result


def map_representatives_to_targets(cluster_entries):
    representative_map = {entry['cluster_id']: entry for entry in cluster_entries if entry['is_representative']}
    pairs_list = []

    for entry in cluster_entries:
        if not entry['is_representative']:
            rep = representative_map[entry['cluster_id']]
            pair_dict = {'representative': rep, 'target': entry}
            pairs_list.append(pair_dict)

    return pairs_list


