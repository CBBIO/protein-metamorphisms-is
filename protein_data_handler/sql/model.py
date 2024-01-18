from sqlalchemy import (Column, Integer, String, Date, ForeignKey, DateTime,
                        func, Float, Boolean)
from sqlalchemy.orm import relationship, declarative_base

Base = declarative_base()


class Protein(Base):
    """
    Represents a protein with its properties and relationships in a database.

    Attributes:
        entry_name (str): Unique entry name for the protein, used as the primary key.
        data_class (str): Classification of the protein's data.
        molecule_type (str): Type of protein molecule.
        sequence_length (int): Length of the protein's amino acid sequence.
        sequence (str): Amino acid sequence of the protein.
        accessions (relationship): Relationship with the 'Accession' class, representing the access codes associated with the protein.
        created_date (Date): Date the protein record was created.
        sequence_update_date (Date): Date of the last update of the protein's sequence.
        annotation_update_date (Date): Date of the last update of the protein's annotation.
        description (str): General description of the protein.
        gene_name (str): Name of the gene associated with this protein.
        organism (str): Organism to which the protein belongs.
        organelle (str): Specific organelle where the protein is located, if applicable.
        organism_classification (str): Taxonomic classification of the organism.
        taxonomy_id (str): Taxonomy identifier for the organism.
        host_organism (str): Host organism of the protein, if applicable.
        host_taxonomy_id (str): Taxonomy identifier for the host organism.
        comments (str): Additional comments about the protein.
        pdb_references (relationship): Relationship with the 'PDBReference' class, containing references to the protein structure database.
        go_terms (relationship): Relationship with the 'GOTerm' class, representing Gene Ontology terms associated with the protein.
        keywords (str): Keywords related to the protein.
        protein_existence (int): Numeric indicator of the protein's existence.
        seqinfo (str): Additional information about the protein's sequence.
        disappeared (Boolean): Indicates whether the protein is no longer present or relevant.
        created_at (DateTime): Date and time the record was created.
        updated_at (DateTime): Date and time of the last update of the record.
    """
    __tablename__ = "proteins"
    entry_name = Column(String, primary_key=True, unique=True, nullable=False)
    data_class = Column(String)
    molecule_type = Column(String)
    sequence = Column(String)
    sequence_length = Column(Integer)
    accessions = relationship(
        "Accession", back_populates="protein"
    )
    created_date = Column(Date)
    sequence_update_date = Column(Date)
    annotation_update_date = Column(Date)
    description = Column(String)
    gene_name = Column(String)
    organism = Column(String)
    organelle = Column(String)
    organism_classification = Column(String)
    taxonomy_id = Column(String)
    host_organism = Column(String)
    host_taxonomy_id = Column(String)
    comments = Column(String)

    pdb_references = relationship("PDBReference", back_populates="protein")
    go_terms = relationship("GOTerm", back_populates="protein")

    keywords = Column(String)
    protein_existence = Column(Integer)
    seqinfo = Column(String)
    disappeared = Column(Boolean)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())


class Accession(Base):
    """
    Represents a unique access code for a protein in a database.

    Attributes:
        id (int): Unique identifier for the access.
        accession_code (str): Unique access code for the protein.
        primary (Boolean): Indicates if the accession code is the primary identifier for the protein.
        protein_entry_name (str): Entry name of the associated protein.
        protein (relationship): Relationship with the 'Protein' class.
        disappeared (Boolean): Indicates whether the access code is no longer in use.
        created_at (DateTime): Date and time the record was created.
        updated_at (DateTime): Date and time of the last update of the record.
    """
    __tablename__ = "accessions"
    id = Column(Integer, primary_key=True)
    accession_code = Column(String, unique=True, nullable=False)
    protein_entry_name = Column(String, ForeignKey("proteins.entry_name"))
    protein = relationship("Protein", back_populates="accessions")
    disappeared = Column(Boolean)
    primary = Column(Boolean)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())


class PDBReference(Base):
    """
    Represents a reference to a structure in the Protein Data Bank (PDB).

    This class is used to store and manage information about protein structures as recorded in the PDB,
    including their relationship to UniProt entries and specific chain details.

    Attributes:
        id (int): Unique identifier for the PDB reference within the database.
        pdb_id (str): Unique identifier of the protein structure in the PDB.
        protein_entry_name (str): Entry name of the associated protein in UniProt.
        protein (relationship): Relationship to the 'Protein' class, linking to the corresponding UniProt entry.
        method (str): Method used to determine the protein structure (e.g., X-ray crystallography, NMR).
        resolution (Float): Resolution of the protein structure, typically in Ångströms (Å).
        uniprot_chains (relationship): Relationship to the 'UniprotChains' class, representing the individual chains in the protein structure as per UniProt.
        pdb_chains (relationship): Relationship to the 'PDBChains' class, detailing the chains in the protein structure as recorded in PDB.
        created_at (DateTime): Timestamp of when the PDB reference record was created in the database.
        updated_at (DateTime): Timestamp of the last update to the PDB reference record.
    """
    __tablename__ = "pdb_references"
    id = Column(Integer, primary_key=True)
    pdb_id = Column(String, nullable=False)
    protein_entry_name = Column(String, ForeignKey("proteins.entry_name"))
    protein = relationship("Protein", back_populates="pdb_references")

    # Method used for the determination of the structure
    method = Column(String)
    resolution = Column(Float)  # Resolution of the structure
    uniprot_chains = relationship("UniprotChains", back_populates="pdb_reference")
    pdb_chains = relationship("PDBChains", back_populates="pdb_reference")
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())


class UniprotChains(Base):
    """
    Represents an individual chain within a protein structure in the Protein Data Bank (PDB) database.

    This class is used to store and manage information about specific chains that are part of a protein structure as
    recorded in the PDB. It includes details about the chain's sequence and its position within the overall protein
    structure.

    Attributes:
        id (int): Unique identifier for each chain in the database.
        pdb_reference_id (int): Foreign key referencing the unique identifier of the protein structure in the PDB
            database to which this chain belongs.
        chain (str): Identifier of the chain within the protein structure, such as 'A', 'B', etc.
        sequence (str): Amino acid sequence of the chain.
        seq_start (int): Starting position of the chain's sequence in the protein structure.
        seq_end (int): Ending position of the chain's sequence in the protein structure.
        pdb_reference (relationship): Relationship with the 'PDBReference' entity, representing the complete protein
            structure to which this chain belongs.
        created_at (DateTime): Date and time when the chain record was created.
        updated_at (DateTime): Date and time of the last update to the chain record.

    The relationship with 'PDBReference' allows each chain to be associated with its corresponding protein structure in
        the PDB database.
    """

    __tablename__ = "uniprot_chains"
    id = Column(Integer, primary_key=True)
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'))
    chain = Column(String)
    sequence = Column(String)
    seq_start = Column(Integer)
    seq_end = Column(Integer)
    pdb_reference = relationship("PDBReference", back_populates="uniprot_chains")
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    def insert_sequence(self, full_sequence):
        """
        Inserta la secuencia en la cadena basándose en seq_start y seq_end.

        Args:
            full_sequence (str): La secuencia completa de la proteína.
        """
        if self.seq_start is not None and self.seq_end is not None:
            self.sequence = full_sequence[self.seq_start - 1:self.seq_end]
        else:
            self.sequence = None


class Cluster(Base):
    """
    Represents a cluster of protein chains, where each cluster is formed by chains with significant similarity,
    determined using the cd-hit tool.

    This class is instrumental in grouping protein chains that are highly similar to each other, aiding in the
    identification of common structures and functions.

    Attributes:
        id (int): Unique identifier for each cluster.
        pdb_chain_id (int): Foreign key referencing the 'PDBChains' entity. It is used to identify the specific protein
        chain in the PDB database associated with this cluster.
        cluster_id (int): Identifier of the cluster, typically a unique string representing this specific group of
            protein chains.
        is_representative (Boolean): Indicates whether the cluster is representative of a larger set of similar chains.
            'True' for yes, 'False' for no.
        sequence_length (int): Average length of the sequences of the chains in the cluster.
        identity (Float): Value representing the average sequence identity within the cluster, usually a percentage
            indicating how similar the chains are within the group.

    The relationship with 'PDBChains' allows each cluster to be connected to its specific chain in the PDB database,
    providing a direct link to detailed structural information.
    """
    __tablename__ = 'clusters'

    id = Column(Integer, primary_key=True)
    pdb_chain_id = Column(Integer, ForeignKey('pdb_chains.id'))
    cluster_id = Column(Integer)
    is_representative = Column(Boolean)
    sequence_length = Column(Integer)
    identity = Column(Float)


class PDBChains(Base):
    """
    Represents an individual chain within a protein structure in the Protein Data Bank (PDB) database.

    Each `PDBChains` object corresponds to a specific polypeptide chain in a protein structure as recorded in the PDB.
    This class is crucial for detailed management of protein structures at the chain level.

    Attributes:
        id (int): Unique identifier for each chain within the database.
        chains (String): Identifier of the chain within the protein structure in PDB. This field is part of the
            composite primary key.
        sequence (String): Amino acid sequence of the protein chain. This field is mandatory and represents the linear
            sequence of amino acids in the chain.
        pdb_reference_id (Integer): Foreign key referencing the unique identifier of the protein structure in the PDB
            database. This field is part of the composite primary key and establishes a direct relationship with the
            `PDBReference` entity.
        pdb_reference (relationship): Relationship with the `PDBReference` entity, providing direct access to detailed
            information about the complete protein structure to which this chain belongs.

    The composite primary key structure of `chains` and `pdb_reference_id` ensures that each instance of `PDBChains` is
    unique and clearly linked to a specific structure in the PDB.
    """
    __tablename__ = 'pdb_chains'
    id = Column(Integer, primary_key=True)
    chains = Column(String)
    sequence = Column(String, nullable=False)
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'))

    pdb_reference = relationship("PDBReference", back_populates="pdb_chains")


class CEAlignResults(Base):
    """
    Represents the results of a CE (Combinatorial Extension) alignment process for protein structures.

    This class stores the results of structural alignment computations, typically involving methods like CEAlign. It is useful for analyzing and comparing the structural alignment of protein clusters.

    Attributes:
        id (int): Unique identifier for each CEAlign result entry in the database.
        cluster_entry_id (int): Foreign key referencing the 'Cluster' entity. It is used to identify the specific protein cluster associated with this alignment result.
        rms (Float): Root Mean Square Deviation (RMSD) value resulting from the CE alignment. RMSD is a measure of the average distance between the atoms (usually the backbone atoms) of superimposed proteins.

    The relationship with 'Cluster' allows each CEAlign result to be directly associated with a specific protein cluster, providing insights into the structural similarity within the cluster.
    """
    __tablename__ = 'ce_align_results'
    id = Column(Integer, primary_key=True)
    cluster_entry_id = Column(Integer, ForeignKey('clusters.id'))
    rms = Column(Float)


class GOTerm(Base):
    """
    Represents a Gene Ontology (GO) term associated with a protein.

    This class is used to store and manage information about the functional annotation of proteins as defined by the Gene Ontology Consortium. Each GO term provides a standardized description of a protein's molecular function, biological process, or cellular component.

    Attributes:
        id (int): Unique identifier for the GO term within the database.
        go_id (str): Unique identifier of the GO term in the Gene Ontology system.
        protein_entry_name (str): Entry name of the associated protein in UniProt.
        protein (relationship): Relationship with the 'Protein' class, linking the GO term to its corresponding protein.
        category (str): Category of the GO term, indicating whether it describes a molecular function, biological process, or cellular component.
        description (str): Detailed description of the GO term, explaining the function, process, or component it represents.

    The relationship with the 'Protein' class allows for the association of functional, process, or component annotations with specific proteins.
    """
    __tablename__ = "go_terms"
    id = Column(Integer, primary_key=True)
    go_id = Column(String, nullable=False)
    protein_entry_name = Column(String, ForeignKey("proteins.entry_name"))
    protein = relationship("Protein", back_populates="go_terms")
    category = Column(String)
    description = Column(String)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())
