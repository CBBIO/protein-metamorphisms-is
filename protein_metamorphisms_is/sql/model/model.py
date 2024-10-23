from protein_metamorphisms_is.sql.model.core.base import Base








class Chains(Base):
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


class StructureEmbedding(Base):
    __tablename__ = 'structure_embeddings'
    id = Column(Integer, primary_key=True)
    model_id = Column(Integer, ForeignKey('models.id'), nullable=False)
    embedding = Column(String)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    model = relationship("Model")
    subcluster_entries = relationship("SubclusterEntry", back_populates="structure_embedding")


class GOAnnotation(Base):
    __tablename__ = 'go_annotation'

    # Primary key
    id = Column(Integer, primary_key=True, autoincrement=True)

    # Foreign keys
    protein_entry_name = Column(String, ForeignKey('proteins.entry_name'), nullable=True)
    go_id = Column(String, ForeignKey('go_terms.go_id'))

    is_transferred = Column(Boolean, default=False)
    embedding_type_id = Column(Integer, ForeignKey('sequence_embedding_types.id'), nullable=True)
    distance = Column(Float, nullable=True)
    source_cluster_id = Column(Integer, ForeignKey('clusters.id'), nullable=True)
    target_cluster_id = Column(Integer, ForeignKey('clusters.id'), nullable=True)

    # Relationships
    protein = relationship("Protein", back_populates="go_term_annotations")
    go_term = relationship("GOTerm")
    source_cluster = relationship("Cluster", foreign_keys=source_cluster_id, back_populates="source_annotations")
    target_cluster = relationship("Cluster", foreign_keys=target_cluster_id, back_populates="target_annotations")
    embedding_type = relationship("SequenceEmbeddingType")

    def __repr__(self):
        return f"<GOAnnotation(id={self.id}, protein={self.protein_entry_name}, go_id={self.go_id}, embedding_type={self.embedding_type_id})>"


class Cluster(Base):
    __tablename__ = 'clusters'
    id = Column(Integer, primary_key=True)
    created_at = Column(DateTime, default=datetime.now)

    # Relación con ClusterEntries
    entries = relationship("ClusterEntry", back_populates="cluster")

    # Relación con GOAnnotation (cluster de origen y destino)
    source_annotations = relationship("GOAnnotation", foreign_keys=GOAnnotation.source_cluster_id,
                                      back_populates="source_cluster")
    target_annotations = relationship("GOAnnotation", foreign_keys=GOAnnotation.target_cluster_id,
                                      back_populates="target_cluster")

    # Adding the missing subclusters relationship
    subclusters = relationship("Subcluster", back_populates="cluster")

    def __repr__(self):
        return f"<Cluster(id={self.id})>"


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


class GOTerm(Base):
    """
    Represents a Gene Ontology (GO) term associated with a protein.
    """
    __tablename__ = "go_terms"

    go_id = Column(String, primary_key=True, nullable=False)
    category = Column(String)
    description = Column(String)
    evidences = Column(String, nullable=True)

    # Relationship with GOPairs
    pairs = relationship("GOPairsTerms", back_populates="go_term")



class GOPairs(Base):
    """
    Represents a pair of GO terms for functional annotation analysis.
    """
    __tablename__ = 'go_pairs'

    # Definición de clave primaria simple
    id = Column(Integer, primary_key=True)
    protein_go_pairs = relationship("ProteinGoPair", back_populates="go_pair")
    # Relación con los resultados de análisis de pares de GO
    results = relationship("GOResultsPairwise", back_populates="pair")


class ProteinGoPair(Base):
    __tablename__ = 'protein_go_pair'

    # Primary key
    id = Column(Integer, primary_key=True, autoincrement=True)

    # Foreign key to GO pairs
    go_pair_id = Column(Integer, ForeignKey('go_pairs.id'), nullable=False)

    # Foreign key to protein
    protein_entry_name = Column(String, ForeignKey('proteins.entry_name'), nullable=False)

    # Relationships
    go_pair = relationship("GOPairs", back_populates="protein_go_pairs")
    protein = relationship("Protein")

    # Restricción de unicidad para evitar duplicados
    __table_args__ = (UniqueConstraint('go_pair_id', 'protein_entry_name', name='_go_pair_protein_uc'),)


class GOPairsTerms(Base):
    """
    Represents the terms associated with a pair of GO terms.
    This class uses a composite primary key combining 'id' and 'go_term_id'.
    """
    __tablename__ = 'go_pairs_terms'

    # Definición de clave primaria compuesta
    id = Column(Integer, ForeignKey('go_pairs.id'), primary_key=True)  # Referencia a la tabla go_pairs
    go_term_id = Column(String, ForeignKey('go_terms.go_id'), primary_key=True)  # Clave foránea como parte de la clave primaria

    # Relación con la clase GOTerm
    go_term = relationship("GOTerm", back_populates="pairs")

    def __repr__(self):
        return f"<GOPairsTerms(id={self.id}, go_term_id={self.go_term_id})>"



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


# class GOPairAnnotations(Base):
#     __tablename__ = 'go_pair_annotations'
#
#     # Primary key
#     id = Column(Integer, primary_key=True, autoincrement=True)
#
#     # Foreign key to GOPairs
#     go_pair_id = Column(Integer, ForeignKey('go_pairs.id'), nullable=False)
#
#     # Foreign keys to two GOAnnotation entries
#     annotation_1_id = Column(Integer, ForeignKey('go_annotation.id'), nullable=False)
#     annotation_2_id = Column(Integer, ForeignKey('go_annotation.id'), nullable=False)
#
#     # Relationships
#     go_pair = relationship("GOPairs", back_populates="annotations")
#     annotation_1 = relationship("GOAnnotation", foreign_keys=[annotation_1_id])
#     annotation_2 = relationship("GOAnnotation", foreign_keys=[annotation_2_id])
