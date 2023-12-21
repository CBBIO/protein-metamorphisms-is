import os

from Bio.PDB import PDBList, PDBParser, CEAligner
import warnings

# Suppress specific PDBConstructionWarning
import warnings
from Bio import PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Suppress specific PDBConstructionWarning
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


from protein_data_handler.models.uniprot import PDBReference, PDBChain, Cluster


class CEAlignHandler:
    def __init__(self, session, download_dir='pdb_files'):
        self.session = session
        self.targets = []
        self.download_dir = download_dir
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)

    def download_structure(self, pdb_id):
        pdbl = PDBList(server='https://files.wwpdb.org/')
        file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=self.download_dir, file_format='pdb')
        return file_path

    def get_targets(self):
        clusters = self.session.query(Cluster).filter_by(is_representative=True).distinct().all()
        self.targets = clusters

    def all_to_all_align(self, identity_threshold):
        parser = PDBParser()
        aligner = CEAligner()

        for target in self.targets:
            representative_chain = self.session.query(PDBChain).filter_by(id=target.pdb_chain_id).one()
            representative_pdb = self.session.query(PDBReference).filter_by(
                id=representative_chain.pdb_reference_id).one()

            # Descargar y parsear la estructura representativa
            rep_file_path = self.download_structure(representative_pdb.pdb_id)
            structure1 = parser.get_structure("rep_structure", rep_file_path)

            clusters_non_representative = self.session.query(Cluster).filter_by(
                is_representative=False, cluster_id=target.cluster_id
            ).filter(Cluster.identity > identity_threshold).all()

            aligner.set_reference(structure1)

            for cluster in clusters_non_representative:
                chain = self.session.query(PDBChain).filter_by(id=cluster.pdb_chain_id).one()
                pdb = self.session.query(PDBReference).filter_by(id=chain.pdb_reference_id).one()

                # Descargar y parsear la estructura no representativa
                non_rep_file_path = self.download_structure(pdb.pdb_id)
                structure2 = parser.get_structure("non_rep_structure", non_rep_file_path)

                for chain_x in structure2.get_chains():
                    print('num atoms',len(list(chain_x.get_atoms())))

                print('-----')
                print(representative_chain.chain, representative_pdb.pdb_id)
                print(chain.chain, pdb.pdb_id)

                # Realizar el alineamiento CE
                # Aquí puedes ajustar los parámetros según tus necesidades
                print(len(chain.sequence),len(representative_chain.sequence))

                print('aligning')
                aligner.align(structure2)
                print(aligner.rms)