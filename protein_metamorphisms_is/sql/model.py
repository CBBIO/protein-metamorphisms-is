from datetime import datetime

from pgvector.sqlalchemy import Vector
from sqlalchemy import (Column, Integer, String, Date, ForeignKey, DateTime,
                        func, Float, Boolean, Index, UniqueConstraint)
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import relationship, declarative_base, mapped_column

Base = declarative_base()


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
    tag = Column(String)
    primary = Column(Boolean)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())


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
        sequence_id (int): Identifier for the associated sequence in the 'Sequence' table.
        sequence (relationship): A link to the 'Sequence' class, indicating the associated sequence.
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
        go_term_associations (relationship): A connection to the 'ProteinGOTermAssociations' class, indicating Gene Ontology terms associated with the protein.
        go_per_protein_semantic_distances (relationship): A connection to the 'GOPerProteinSemanticDistance' class, indicating semantic distances for GO terms per protein.
        keywords (str): Descriptive keywords related to the protein, aiding in categorization and search.
        protein_existence (int): A numerical code indicating the evidence level for the protein's existence.
        seqinfo (str): Supplementary information about the protein's sequence.
        disappeared (Boolean): Flag indicating whether the protein is obsolete or no longer relevant.
        created_at (DateTime): Timestamp of when the record was initially created.
        updated_at (DateTime): Timestamp of the most recent update to the record.
        structure_id (int): Identifier for the associated structure file in the 'Structure' table.
        structure (relationship): A link to the 'Structure' class, indicating the associated structure file.
    """
    __tablename__ = "proteins"
    entry_name = Column(String, primary_key=True, unique=True, nullable=False)
    data_class = Column(String)
    molecule_type = Column(String)
    sequence_id = Column(Integer, ForeignKey('sequences.id'))
    sequence = relationship("Sequence", uselist=False)
    accessions = relationship("Accession", back_populates="protein")
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
    go_term_associations = relationship("ProteinGOTermAssociations", back_populates="protein")
    go_per_protein_semantic_distances = relationship("GOPerProteinSemanticDistance", back_populates="protein")
    keywords = Column(String)
    protein_existence = Column(Integer)
    seqinfo = Column(String)
    disappeared = Column(Boolean)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())
    structure_id = Column(Integer, ForeignKey('structures.id'), nullable=True)
    structure = relationship("Structure")


class Sequence(Base):
    __tablename__ = 'sequences'
    id = Column(Integer, primary_key=True)
    sequence = Column(String, nullable=False)
    sequence_hash = Column(String, index=True, unique=True)

    # Adding a default value to automatically compute the hash when a sequence is added
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sequence_hash = func.md5(self.sequence)


Index('idx_sequence_hash', Sequence.sequence_hash)


class Structure(Base):
    """
    Represents a structure file (.cif) and its associated metadata.

    Attributes:
        id (int): Unique identifier for the structure.
        file_path (str): Relative path to the structure file (.cif).
        created_at (DateTime): Timestamp of when the record was initially created.
        updated_at (DateTime): Timestamp of the most recent update to the record.
    """
    __tablename__ = 'structures'
    id = Column(Integer, primary_key=True)
    file_path = Column(String, nullable=False)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())
    models = relationship("Model", back_populates="structure")  # Relación con Model


class Model(Base):
    """
    Represents different models within a structure, specifically for NMR structures that may include multiple models. Each model
    is identified by a model_id that is typically assigned during the experimental crystallographic process.

    Attributes:
        id (int): Unique identifier for each model within a structure.
        model_id (str): Identifier for the model as assigned in crystallography.
        structure_id (int): Foreign key to the structure this model belongs to.
        structure (relationship): A relationship back to the 'Structure' class.
    """
    __tablename__ = 'models'
    id = Column(Integer, primary_key=True)
    model_id = Column(Integer, nullable=False)
    file_path = Column(String, nullable=False)
    structure_id = Column(Integer, ForeignKey('structures.id'), nullable=False)
    structure = relationship("Structure", back_populates="models")


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
    # sequence_id = Column(Integer, ForeignKey('sequences.id'))
    # sequence = relationship("Sequence", uselist=False)  # Single sequence per protein
    resolution = Column(Float)
    pdb_chains = relationship("PDBChains", back_populates="pdb_reference")
    structure_id = Column(Integer, ForeignKey('structures.id'), nullable=True)
    structure = relationship("Structure")
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
    sequence_id = Column(Integer, ForeignKey('sequences.id'))
    sequence = relationship("Sequence", uselist=False)
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'))
    structure_id = Column(Integer, ForeignKey('structures.id'), nullable=False)  # Añadido para vincular con Structure
    structure = relationship("Structure")
    pdb_reference = relationship("PDBReference", back_populates="pdb_chains")


class SequenceEmbedding(Base):
    __tablename__ = 'sequence_embeddings'
    id = Column(Integer, primary_key=True)
    sequence_id = Column(Integer, ForeignKey('sequences.id'), nullable=False)
    embedding_type_id = Column(Integer, ForeignKey('sequence_embedding_types.id'))
    embedding = mapped_column(Vector())
    shape = Column(ARRAY(Integer))  # Almacena las dimensiones del embedding como un array de enteros
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    sequence = relationship("Sequence")
    embedding_type = relationship("SequenceEmbeddingType")

    # Relación con las predicciones de GO
    go_predictions = relationship("SequenceEmbeddingGOAnnotationTransfer", back_populates="embedding")




class StructureEmbedding(Base):
    __tablename__ = 'structure_embeddings'
    id = Column(Integer, primary_key=True)
    model_id = Column(Integer, ForeignKey('models.id'), nullable=False)
    embedding = Column(String)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    model = relationship("Model")
    subcluster_entries = relationship("SubclusterEntry", back_populates="structure_embedding")


class Cluster(Base):
    __tablename__ = 'clusters'
    id = Column(Integer, primary_key=True)
    created_at = Column(DateTime, default=datetime.now)
    # Relación con ClusterEntries
    entries = relationship("ClusterEntry", back_populates="cluster")
    # Relación con Subcluster
    subclusters = relationship("Subcluster", back_populates="cluster")


class ClusterEntry(Base):
    __tablename__ = 'cluster_entries'
    id = Column(Integer, primary_key=True)
    cluster_id = Column(Integer, ForeignKey('clusters.id'))
    sequence_id = Column(Integer, ForeignKey('sequences.id'), nullable=False)
    is_representative = Column(Boolean)
    sequence_length = Column(Integer)
    identity = Column(Float)
    created_at = Column(DateTime, default=datetime.now)

    # Relaciones con Cluster y PDBChains
    cluster = relationship("Cluster", back_populates="entries")
    sequence = relationship("Sequence")


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


class SequenceEmbeddingType(Base):
    """
    Represents a type of protein sequence analysis embedding.

    This class is designed to manage different embedding techniques used in protein sequence analysis, offering a structured way to categorize and store information about various embedding methods such as ESM and Prot-T5.

    Attributes:
        id (Integer): Unique identifier for each embedding type.
        name (String): Unique name of the embedding type.
        description (String): Detailed description of the embedding technique.
        task_name (String): Name of the specific task associated with this embedding type, if applicable.
    """
    __tablename__ = 'sequence_embedding_types'
    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False, unique=True)
    description = Column(String)
    task_name = Column(String)
    model_name = Column(String)

    seq_embeddings = relationship("SequenceEmbedding", back_populates="embedding_type")


class StructureEmbeddingType(Base):
    """
    Represents a type of protein structure analysis embedding.

    This class is designed to manage different embedding techniques used in protein structure analysis, offering a structured way to categorize and store information about various embedding methods such as 3D-CNN and GNN.

    Attributes:
        id (Integer): Unique identifier for each embedding type.
        name (String): Unique name of the embedding type.
        description (String): Detailed description of the embedding technique.
        task_name (String): Name of the specific task associated with this embedding type, if applicable.
    """
    __tablename__ = 'structure_embedding_types'
    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False, unique=True)
    description = Column(String)
    task_name = Column(String)
    model_name = Column(String)


class StructuralAlignmentResults(Base):
    """
    Stores results from structural alignment tasks, with ordered pairs of SubclusterEntry.

    Attributes:
        id (Integer): Unique identifier for each result set.
        cluster_id (Integer, ForeignKey): Reference to the parent cluster.
        subcluster_1_id (Integer, ForeignKey): Reference to the first SubclusterEntry in the ordered pair.
        subcluster_2_id (Integer, ForeignKey): Reference to the second SubclusterEntry in the ordered pair.
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
    cluster_id = Column(Integer, ForeignKey('clusters.id'))
    subcluster_1_id = Column(Integer, ForeignKey('subcluster_entries.id'), nullable=False)
    subcluster_2_id = Column(Integer, ForeignKey('subcluster_entries.id'), nullable=False)

    # Alignment results
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


class ProteinGOTermAssociations(Base):
    __tablename__ = 'protein_go_term_association'
    protein_entry_name = Column(String, ForeignKey('proteins.entry_name'), primary_key=True)
    go_id = Column(String, ForeignKey('go_terms.go_id'), primary_key=True)
    protein = relationship(
        "Protein",
        back_populates="go_term_associations",
    )
    go_term = relationship(
        "GOTerm"
    )


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
    go_id = Column(String, primary_key=True, nullable=False)
    category = Column(String)
    description = Column(String)
    evidences = Column(String, nullable=True)
    sequence_predictions = relationship("SequenceEmbeddingGOAnnotationTransfer", back_populates="go_term")


class GOPairs(Base):
    """
    Represents pairs of Gene Ontology (GO) terms associated with a specific protein.

    Attributes:
        id (int): Unique identifier for the pair.
        go_term_1_id (str): GO term ID for the first GO term in the pair.
        go_term_2_id (str): GO term ID for the second GO term in the pair.
    """
    __tablename__ = 'go_pairs'
    id = Column(Integer, primary_key=True)
    go_term_1_id = Column(String, ForeignKey('go_terms.go_id'))
    go_term_2_id = Column(String, ForeignKey('go_terms.go_id'))

    # Relaciones
    go_term_1 = relationship("GOTerm", foreign_keys=[go_term_1_id])
    go_term_2 = relationship("GOTerm", foreign_keys=[go_term_2_id])

    results = relationship("GOResultsPairwise", back_populates="pair")


class GOResultsPairwise(Base):
    """
    Represents the results of analysis on pairs of GO terms, storing distances and other metrics.

    Attributes:
        id (int): Unique identifier for the result.
        pair_id (int): Reference to the associated `GOPairs` entry.
        information_content_1 (float): Information content of the first GO term.
        information_content_2 (float): Information content of the second GO term.
        resnik_distance (float): Resnik distance between the two GO terms.
        minimum_branch_length (float): Minimum branch length between the two GO terms.
    """
    __tablename__ = 'go_results_pairwise'
    id = Column(Integer, primary_key=True)
    pair_id = Column(Integer, ForeignKey('go_pairs.id'))
    information_content_1 = Column(Float)
    information_content_2 = Column(Float)
    resnik_distance = Column(Float)
    minimum_branch_length = Column(Float)

    # Relación con los pares de GO
    pair = relationship("GOPairs", back_populates="results")


class SequenceEmbeddingGOAnnotationTransfer(Base):
    __tablename__ = 'sequence__embedding_go_annotation_transfer'
    id = Column(Integer, primary_key=True)
    embedding_id = Column(Integer, ForeignKey('sequence_embeddings.id'))
    ref_protein_entry_name = Column(String, ForeignKey('proteins.entry_name'))
    go_id = Column(String, ForeignKey('go_terms.go_id'), nullable=False)
    embedding_type_id = Column(Integer, ForeignKey('sequence_embedding_types.id'), nullable=False)
    prediction_method_id = Column(Integer, ForeignKey('prediction_methods.id'))
    k = Column(Integer)

    embedding = relationship("SequenceEmbedding", back_populates="go_predictions")
    go_term = relationship("GOTerm", back_populates="sequence_predictions")
    prediction_method = relationship("PredictionMethod")



class PredictionMethod(Base):
    __tablename__ = 'prediction_methods'
    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False, unique=True)
    description = Column(String)


class GOPerProteinSemanticDistance(Base):
    __tablename__ = 'go_per_protein_semantic_distances'
    id = Column(Integer, primary_key=True)
    protein_entry_name = Column(String, ForeignKey('proteins.entry_name'))
    embedding_type_id = Column(Integer, ForeignKey('sequence_embedding_types.id'))
    prediction_method_id = Column(Integer, ForeignKey('prediction_methods.id'))  # Vinculación con PredictionMethod
    group_distance = Column(Float, nullable=False)

    protein = relationship("Protein", back_populates="go_per_protein_semantic_distances")
    embedding_type = relationship("SequenceEmbeddingType")
    prediction_method = relationship("PredictionMethod")


class Subcluster(Base):
    __tablename__ = 'subclusters'
    id = Column(Integer, primary_key=True)
    cluster_id = Column(Integer, ForeignKey('clusters.id'))
    description = Column(String, nullable=True)
    created_at = Column(DateTime, default=datetime.now)

    # Relaciones
    cluster = relationship("Cluster", back_populates="subclusters")
    entries = relationship("SubclusterEntry", back_populates="subcluster")


class SubclusterEntry(Base):
    __tablename__ = 'subcluster_entries'
    id = Column(Integer, primary_key=True)
    subcluster_id = Column(Integer, ForeignKey('subclusters.id'))
    structure_embedding_id = Column(Integer, ForeignKey('structure_embeddings.id'))
    is_representative = Column(Boolean, default=False)
    created_at = Column(DateTime, default=datetime.now)

    # Relaciones
    subcluster = relationship("Subcluster", back_populates="entries")
    structure_embedding = relationship("StructureEmbedding", back_populates="subcluster_entries")
