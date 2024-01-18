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

    def start(self):
        """
        Begins the structural alignment process. This method orchestrates the entire alignment workflow,
        starting from loading the cluster representatives, performing CE alignment, and handling any exceptions
        that occur during the process.
        """
        try:
            cluster_representatives = self.load_clusters()
            self.ce_align(cluster_representatives)

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
        cluster_representatives = self.session.query(Cluster).filter_by(is_representative=True).all()
        return cluster_representatives

    def ce_align(self, cluster_representatives):
        """
        Performs the CE alignment on the cluster representatives. This method utilizes concurrent processing
        to align multiple protein structures in parallel, significantly speeding up the computation time.

        Args:
            cluster_representatives (list): List of cluster representative objects to align.
        """
        num_workers = self.conf.get('num_workers', 4)
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(self.align_task, rep) for rep in cluster_representatives]

            for future in as_completed(futures):
                try:
                    alignment_results = future.result()
                except Exception as e:
                    self.logger.error(f"Error in alignment task: {e}")

    def align_task(self, representative):
        """
        Task to align a single representative against all targets in its cluster. This function is where
        the actual structural alignment takes place using the CE algorithm. It includes loading the protein
        structures from PDB files, aligning them, and computing the RMS (Root Mean Square) deviation,
        which quantifies the similarity between the structures.

        Args:
            representative (object): The representative object to be aligned.
        """

        pdb_chains_path = self.conf.get('pdb_chains_path', './chains')
        parser = PDBParser()

        representative_details = self.get_cluster_pdb_chain_details(representative.id)
        representative_name = f"{representative_details[2]}_{representative_details[1]}"
        representative_structure_path = os.path.join(pdb_chains_path, f"{representative_name}.pdb")

        representative_structure = parser.get_structure(representative_name, representative_structure_path)

        targets = self.session.query(Cluster).filter_by(cluster_id=representative.cluster_id,
                                                        is_representative=False).all()

        aligner = CEAligner()
        aligner.set_reference(representative_structure)

        local_session = sessionmaker(bind=self.engine)()

        for target in targets:
            target_details = self.get_cluster_pdb_chain_details(target.id)
            target_name = f"{target_details[2]}_{target_details[1]}"
            target_structure_path = os.path.join(pdb_chains_path, f"{target_name}.pdb")
            target_structure = parser.get_structure(target_name, target_structure_path)

            aligner.align(target_structure)

            rms = aligner.rms

            result = CEAlignResults(cluster_entry_id=target.id, rms=rms)
            local_session.add(result)
        local_session.commit()
        local_session.close()

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

        return result
