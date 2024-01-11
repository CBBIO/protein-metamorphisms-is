import os
from concurrent.futures import ProcessPoolExecutor
from functools import partial

from Bio.PDB import PDBList, MMCIFParser, MMCIFIO, CEAligner
from Bio.PDB.Structure import Structure

from protein_data_handler.sql.model import PDBChains, PDBChain, PDBReference


class ChainExtractor:
    def __init__(self, session, pdb_dir, chain_dir):
        self.session = session
        self.pdb_dir = pdb_dir
        self.chain_dir = chain_dir
        print('extract')

    def decompose(self):
        chains = self.get_chains()
        for chain_group in chains:
            chains_id = chain_group.id
            for chain in chain_group.chains.split('/'):
                pdb_chain = PDBChain(chain=chain, pdb_chains_id=chains_id)
                self.session.add(pdb_chain)
        self.session.commit()

    def get_chains(self):
        return self.session.query(PDBChains).all()

    def create_decomposed_structure(self):
        pdbl = PDBList(server='https://files.wwpdb.org/')
        mmcif_parser = MMCIFParser()

        # Realizar la consulta con join incluyendo PDBReference.pdb_id
        pdb_chains = self.session.query(
            PDBChain, PDBChains.pdb_reference_id, PDBReference.id, PDBReference.pdb_id
        ).join(
            PDBChains, PDBChain.pdb_chains_id == PDBChains.id
        ).join(
            PDBReference, PDBChains.pdb_reference_id == PDBReference.id
        ).all()[:10]

        structures_path = []

        for pdb_chain, pdb_reference_id, pdb_reference_id, pdb_id in pdb_chains:
            print(f"PDB Chain: {pdb_chain.chain}, PDB ID: {pdb_id}")

            # Descargar el archivo PDB
            pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, pdir=self.pdb_dir)
            # Parsear la estructura PDB
            structure = mmcif_parser.get_structure(pdb_id, pdb_file_path)
            new_structure = extract_chain(structure, pdb_chain)

            io = MMCIFIO()
            io.set_structure(new_structure)

            structures_path.append(os.path.join(self.chain_dir, f"{pdb_id}_{pdb_chain.chain}.cif"))

            io.save(os.path.join(self.chain_dir, f"{pdb_id}_{pdb_chain.chain}.cif"))
        reference_structure = mmcif_parser.get_structure('nan', structures_path[-1])
        ce_aligner = CEAligner()
        ce_aligner.set_reference(reference_structure)

        # Paralelizar el cálculo de RMSD
        with ProcessPoolExecutor() as executor:
            # Crear una función parcial con los argumentos constantes
            func = partial(align_structure, reference_structure=reference_structure, ce_aligner=ce_aligner)

            # Utilizar executor.map para aplicar la función a cada estructura
            rmsd_values = executor.map(func, structures_path[:-1])
            for rmsd_value in rmsd_values:
                print(rmsd_value)



def extract_chain(structure, chain_id):
    # Crear una nueva estructura
    new_structure = Structure(structure.id)

    # Añadir el modelo (usualmente hay un solo modelo en estructuras cristalográficas)
    for model in structure:
        new_model = model.copy()
        new_structure.add(new_model)

        # Añadir solo la cadena de interés
        for chain in model:
            if chain.id == chain_id:
                new_chain = chain.copy()
                new_model.add(new_chain)
                break

    return new_structure


def align_structure(self, structure_path, reference_structure, ce_aligner):
    mmcif_parser = MMCIFParser()
    structure = mmcif_parser.get_structure('nan', structure_path)
    ce_aligner.align(structure)
    return ce_aligner.rms