from sqlalchemy import (Column, Integer, String, Date, ForeignKey, DateTime,
                        func, Float, Boolean)
from sqlalchemy.orm import relationship, declarative_base

Base = declarative_base()


class Protein(Base):
    """
    Represents a protein, encapsulating its properties and relationships within a database.

    This class models a protein entity, encompassing various attributes that describe its
    characteristics and relationships to other entities. It serves as a comprehensive record
    for proteins, covering aspects from basic sequence data to more complex annotations and references.

    Attributes:
        entry_name (str): Unique entry name for the protein, serving as the primary key.
        data_class (str): Categorization of the protein's data (e.g., experimental, predicted).
        molecule_type (str): Type of the protein molecule (e.g., enzyme, antibody).
        sequence_length (int): The length of the amino acid sequence of the protein.
        sequence (str): Full amino acid sequence of the protein.
        accessions (relationship): A link to the 'Accession' class, detailing access codes associated with this protein.
        created_date (Date): The date when the protein record was first created.
        sequence_update_date (Date): The date when the protein's sequence was last updated.
        annotation_update_date (Date): The date when the protein's annotation was last updated.
        description (str): A general description or overview of the protein.
        gene_name (str): The name of the gene that encodes this protein.
        organism (str): The organism from which the protein is derived.
        organelle (str): The specific organelle where the protein is localized, if applicable.
        organism_classification (str): Taxonomic classification of the organism (e.g., species, genus).
        taxonomy_id (str): A unique identifier for the organism in taxonomic databases.
        host_organism (str): The host organism for the protein, relevant in cases of viral or symbiotic proteins.
        host_taxonomy_id (str): Taxonomy identifier for the host organism, if applicable.
        comments (str): Additional remarks or notes about the protein.
        pdb_references (relationship): A link to the 'PDBReference' class, providing references to structural data in the PDB.
        go_terms (relationship): A connection to the 'GOTerm' class, indicating Gene Ontology terms associated with the protein.
        keywords (str): Descriptive keywords related to the protein, aiding in categorization and search.
        protein_existence (int): A numerical code indicating the evidence level for the protein's existence.
        seqinfo (str): Supplementary information about the protein's sequence.
        disappeared (Boolean): Flag indicating whether the protein is obsolete or no longer relevant.
        created_at (DateTime): Timestamp of when the record was initially created.
        updated_at (DateTime): Timestamp of the most recent update to the record.

    This class is integral to managing and querying detailed protein data, supporting a wide range of bioinformatics
    and data analysis tasks.
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
    Represents a unique access code for a protein in a database, often used in bioinformatics repositories.

    This class models an accession record, which is essential for tracking and referencing protein data. Each accession
    record provides a unique identifier for a protein and is linked to detailed protein information.

    The `Accession` class plays a crucial role in the organization and retrieval of protein data, acting as a key
    reference point for protein identification and database querying.

    Attributes:
        id (int): A unique identifier for the accession record within the database.
        accession_code (str): The unique access code associated with a specific protein. This code is typically used as
                              a reference in various bioinformatics databases and literature.
        primary (Boolean): A flag indicating whether this accession code is the primary identifier for the associated protein.
                           Primary accession codes are generally the most stable and widely used references.
        protein_entry_name (str): The entry name of the protein associated with this accession code. This serves as a link
                                  to the protein's detailed record.
        protein (relationship): A SQLAlchemy relationship with the 'Protein' class. This relationship provides a direct
                                connection to the protein entity that this accession code represents, allowing for the retrieval
                                of comprehensive protein information.
        disappeared (Boolean): A flag indicating whether the accession code is obsolete or no longer in use. This is important
                               for maintaining the integrity and relevance of the database.
        created_at (DateTime): The date and time when this accession record was first created in the database.
        updated_at (DateTime): The date and time when this accession record was last updated, reflecting any changes or
                               updates to the accession information.
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

    This class is pivotal for storing and managing details about protein structures as cataloged in the PDB. It forms a bridge
    between PDB structures and UniProt entries, enabling comprehensive tracking and analysis of protein structures and their
    corresponding sequences.

    The `PDBReference` class serves as a critical component for integrating structural data with protein sequence and functional information,
    thereby enriching the understanding of protein structures.

    Attributes:
        id (int): A unique identifier for the PDB reference within the database. This serves as the primary key.
        pdb_id (str): The unique identifier of the protein structure in PDB, typically a 4-character alphanumeric code.
        protein_entry_name (str): The entry name of the associated protein in UniProt. This helps link the structure to its
                                  corresponding protein sequence and other relevant data in UniProt.
        protein (relationship): A SQLAlchemy relationship to the 'Protein' class, establishing a connection to the UniProt entry
                                corresponding to this PDB structure.
        method (str): The method used for determining the protein structure, such as X-ray crystallography or NMR spectroscopy.
        resolution (Float): The resolution of the protein structure, measured in Ångströms (Å). A lower number indicates higher resolution.
        uniprot_chains (relationship): A relationship to the 'UniprotChains' class, detailing the individual protein chains in the structure
                                       as defined in UniProt.
        pdb_chains (relationship): A relationship to the 'PDBChains' class, describing the chains in the protein structure as recorded in PDB.
        created_at (DateTime): The timestamp indicating when the PDB reference record was initially created in the database.
        updated_at (DateTime): The timestamp of the most recent update to the PDB reference record. This field is automatically updated
                               on each record modification.
    """
    __tablename__ = "pdb_references"
    id = Column(Integer, primary_key=True)
    pdb_id = Column(String, nullable=False)
    protein_entry_name = Column(String, ForeignKey("proteins.entry_name"))
    protein = relationship("Protein", back_populates="pdb_references")
    method = Column(String)
    resolution = Column(Float)
    pdb_chains = relationship("PDBChains", back_populates="pdb_reference")
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())


class PDBChains(Base):
    """
    Represents an individual polypeptide chain within a protein structure as cataloged in the Protein Data Bank (PDB).

    The `PDBChains` class is instrumental in representing each distinct polypeptide chain encountered in protein structures
    from the PDB. This class enables detailed tracking and management of these chains, facilitating analyses and queries at
    the chain level. By associating each chain with its parent protein structure, the class enhances the database's ability
    to model complex protein structures.

    Attributes:
        id (int): A unique identifier for each polypeptide chain within the database, serving as the primary key.
        chains (String): The specific identifier of the chain as referenced in the protein structure within PDB. This attribute,
                         combined with 'pdb_reference_id', constitutes part of the composite primary key.
        sequence (String): The complete amino acid sequence of the chain. Storing this mandatory attribute allows for in-depth
                           analyses of the chain's molecular structure.
        pdb_reference_id (Integer): A foreign key linking to the unique identifier of the parent protein structure in the PDB.
                                    This attribute forms the other part of the composite primary key and establishes a direct
                                    relationship with the `PDBReference` entity.
        model (Integer): An identifier for the model of the chain, particularly important for structures like NMR that may
                         encompass multiple models.
        pdb_reference (relationship): A SQLAlchemy relationship that connects to the `PDBReference` entity. This relationship
                                      provides access to comprehensive details about the entire protein structure to which this
                                      chain is a part.

    The composite primary key, comprising `chains` and `pdb_reference_id`, ensures that each instance of `PDBChains` is uniquely
    tied to a specific structure in the PDB. This key structure is critical for precise data retrieval and efficient management
    of the database's structural data.
    """
    __tablename__ = 'pdb_chains'
    id = Column(Integer, primary_key=True)
    chains = Column(String)
    sequence = Column(String, nullable=False)
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'))
    model = Column(Integer)

    pdb_reference = relationship("PDBReference", back_populates="pdb_chains")


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
    complexity_level_id = Column(Integer, ForeignKey('structural_complexity_levels.id'))

    complexity_level = relationship("StructuralComplexityLevel", backref="clusters")


class StructuralComplexityLevel(Base):
    """
    Captures the hierarchy of structural forms within proteins, ranging from individual proteins to the partitioning of
    chains through its secondary structure.

    This class provides a foundational abstraction for handling proteins at various levels of structural complexity within
    the development environment. It allows for the execution of operations across different complexity levels, enabling a
    more flexible and nuanced approach to protein data manipulation and analysis. By defining distinct levels of structural
    complexity, it supports targeted queries and operations, enhancing the efficiency and precision of bioinformatics
    workflows.

    Attributes:
        id (Integer): Unique identifier for each complexity level.
        name (String): Descriptive name of the complexity level.
        description (String, optional): More detailed information about the complexity level.
    """
    __tablename__ = 'structural_complexity_levels'
    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    description = Column(String, nullable=True)


class StructuralAlignmentType(Base):
    """
    Provides a framework for aligning protein structures, crucial for understanding
    the functional and evolutionary relationships between proteins. This class
    enables the use of various alignment strategies, supporting a comprehensive
    approach to protein comparison.

    Structural alignment methods integrated within this framework include:

    - CE-align: Identifies optimal alignments based on the Combinatorial Extension
      method, focusing on similar backbone arrangements.
    - US-align: Utilizes an advanced algorithm for measuring structural similarity,
      offering insights into sequence identity and alignment scores.
    - FATCAT: Capable of accommodating protein flexibility during alignment,
      allowing for the detection of functionally important variations.

    By incorporating these methodologies, the class facilitates diverse approaches
    to protein comparison. This enables researchers to gain deeper insights into
    protein functionality and evolution, highlighting the significance of structural
    alignment in the field of bioinformatics.

    Attributes:
        id (Integer): Unique identifier for each alignment type.
        name (String): Name of the alignment type.
        description (String): Detailed description of the alignment method.
        task_name (String): Name of the specific task or process associated with this alignment type.
    """
    __tablename__ = 'structural_alignment_types'
    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    description = Column(String)
    task_name = Column(String)


class StructuralAlignmentQueue(Base):
    """
        Manages a queue of pending structural alignment tasks, overseeing their execution and monitoring.

        Attributes:
            id (Integer): Unique identifier for each queue entry.
            cluster_entry_id (Integer, ForeignKey): Reference to the protein chain cluster being aligned.
            alignment_type_id (Integer, ForeignKey): ID of the structural alignment type to be applied.
            state (Integer): Current state of the task (e.g., pending, processing, completed, error).
            retry_count (Integer): Number of retries attempted for the task.
            error_message (String, optional): Error message if the task fails.
            created_at (DateTime): Timestamp when the queue entry was created.
            updated_at (DateTime): Timestamp when the queue entry was last updated.
    """
    __tablename__ = 'structural_alignment_queue'
    id = Column(Integer, primary_key=True)
    cluster_entry_id = Column(Integer, ForeignKey('clusters.id'), nullable=False)
    alignment_type_id = Column(Integer, ForeignKey('structural_alignment_types.id'), nullable=False)
    state = Column(Integer, default=0, nullable=False)
    retry_count = Column(Integer, default=0, nullable=False)
    error_message = Column(String, nullable=True)
    created_at = Column(DateTime, default=func.now(), nullable=False)
    updated_at = Column(DateTime, default=func.now(), onupdate=func.now())

    alignment_type = relationship("StructuralAlignmentType")


class StructuralAlignmentResults(Base):
    """
    Stores results from structural alignment tasks, providing detailed metrics and scores.

    Attributes:
        id (Integer): Unique identifier for each set of results.
        cluster_entry_id (Integer, ForeignKey): Reference to the cluster of protein chains analyzed.
        ce_rms (Float): Root mean square deviation calculated by CE method.
        tm_rms (Float): Root mean square deviation calculated by US-align.
        tm_seq_id (Float): Sequence identity calculated by US-align.
        tm_score_chain_1 (Float): TM score for the first chain in the US-alignment.
        tm_score_chain_2 (Float): TM score for the second chain in the US-alignment.
        fc_rms (Float): Root mean square deviation calculated by FATCAT.
        fc_identity (Float): Sequence identity calculated by FATCAT.
        fc_similarity (Float): Similarity score calculated by FATCAT.
        fc_score (Float): Overall score calculated by FATCAT.
        fc_align_len (Float): Length of the alignment calculated by FATCAT.
    """
    __tablename__ = 'structural_alignment_results'
    id = Column(Integer, primary_key=True)
    cluster_entry_id = Column(Integer, ForeignKey('clusters.id'))
    ce_rms = Column(Float)
    tm_rms = Column(Float)
    tm_seq_id = Column(Float)
    tm_score_chain_1 = Column(Float)
    tm_score_chain_2 = Column(Float)
    fc_rms = Column(Float)
    fc_identity = Column(Float)
    fc_similarity = Column(Float)
    fc_score = Column(Float)
    fc_align_len = Column(Float)


class GOTerm(Base):
    """
    Represents a Gene Ontology (GO) term associated with a protein.

    This class is used to store and manage information about the functional annotation of proteins as defined by the Gene Ontology Consortium.
    Each GO term provides a standardized description of a protein's molecular function, biological process, or cellular component.

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
