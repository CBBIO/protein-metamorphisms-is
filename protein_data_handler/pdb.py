from Bio.PDB import PDBList


def download_entire_pdb(server, pdb, file_format):
    """
    Descarga el conjunto completo de estructuras de proteínas en formato mmCIF.

    Esta función utiliza `PDBList` de Biopython para descargar todas las
    estructuras de proteínas disponibles en el Protein Data Bank (PDB).
    Las estructuras se descargan en formato mmCIF y se guardan en un
    directorio local.

    :return: None
    """
    pdb_list = PDBList(server=server, pdb=pdb)
    pdb_list.download_entire_pdb(file_format=file_format)
