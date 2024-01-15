import os
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio.PDB import PDBParser, CEAligner
from sqlalchemy.orm import sessionmaker

from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.operations.base.bioinfo_operator import BioinfoOperatorBase
from protein_data_handler.sql.model import PDBChains, Cluster, PDBReference, CEAlignResults


class CEAlign(BioinfoOperatorBase):
    def __init__(self, conf):
        super().__init__(conf, session_required=True)
        self.num_workers = conf.get('num_workers', 4)  # Número de hilos de ejecución en paralelo

    def start(self):
        try:
            cluster_representatives = self.load_clusters()
            self.ce_align(cluster_representatives)

        except Exception as e:
            self.logger.error(f"Error during structural alignment process: {e}")
            raise

    def load_clusters(self):
        cluster_representatives = self.session.query(Cluster).filter_by(is_representative=True).all()
        return cluster_representatives

    def ce_align(self, cluster_representatives):
        with ThreadPoolExecutor(max_workers=self.num_workers) as executor:
            futures = [executor.submit(self.align_task, rep) for rep in cluster_representatives]

            for future in as_completed(futures):
                try:
                    alignment_results = future.result()
                except Exception as e:
                    self.logger.error(f"Error in alignment task: {e}")

    def align_task(self, representative):

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

    def save_alignment_results(self, results):
        for cluster_id, rms in results:
            result = CEAlignResults(cluster_entry_id=cluster_id, rms=rms)
            self.session.add(result)

        try:
            self.session.commit()
        except Exception as e:
            self.logger.error(f"Error al guardar los resultados de alineación: {e}")
            self.session.rollback()


