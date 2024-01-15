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
    Represents a reference to the Protein Data Bank (PDB) database.

    Attributes:
        id (int): Unique identifier for the PDB reference.
        pdb_id (str): Unique identifier in the PDB database.
        protein_entry_name (str): Entry name of the associated protein.
        protein (relationship): Relationship with the 'Protein' class.
        method (str): Method used for determining the protein structure.
        resolution (Float): Resolution of the protein structure in the PDB database.
        uniprot_chains (relationship): Relationship with the 'UniprotChains' class, representing Uniprot chains in the
         protein structure.
        pdb_chains (relationship): Relationship with the 'PDBChains' class, representing the structure chains in the
         PDB database.
        uniprot_pdb_alignments (relationship): Relationship with the 'UniProtPDBAlignment' class, representing alignments
         between Uniprot and PDB sequences.
        created_at (DateTime): Date and time the record was created.
        updated_at (DateTime): Date and time of the last update of the record.
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
    uniprot_pdb_alignments = relationship("UniProtPDBAlignment", back_populates="pdb_reference")
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())



class UniprotChains(Base):
    """
    Representa una cadena individual dentro de una estructura de proteína
    en la base de datos PDB.

    Esta clase se utiliza para almacenar información sobre cadenas
    específicas que forman parte de una estructura de proteína, como
    se registra en la base de datos de estructuras de proteínas (PDB).

    Attributes:
        id (int): Identificador único para cada cadena en la base de datos.
        pdb_reference_id (int): Clave foránea que referencia al
            identificador único de la estructura de proteína en la base de
            datos PDB a la que pertenece esta cadena.
        chain (str): Identificador de la cadena dentro de la estructura
            de proteína.
            Por ejemplo, 'A', 'B', etc.
        seq_start (int): Posición inicial de la secuencia de la cadena en la
            estructura de proteína.
        seq_end (int): Posición final de la secuencia de la cadena en
            la estructura de proteína.
        pdb_reference (relationship): Relación con la entidad 'PDBReference'
            que representa la estructura completa de la proteína a la que
            pertenece esta cadena.
        created_at (DateTime): Fecha y hora de creación del registro de
            la cadena.
        updated_at (DateTime): Fecha y hora de la última actualización del
            registro de la cadena.

    La relación con 'PDBReference' permite asociar cada cadena con su
    estructura de proteína correspondiente en la base de datos PDB.
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
        Representa un cluster de cadenas de proteínas, donde cada cluster está formado por cadenas con similitud
         significativa, determinado mediante el uso de la herramienta cd-hit.

        Esta clase es útil para agrupar cadenas de proteínas que son muy similares entre sí, lo que ayuda en la
         identificación de estructuras y funciones comunes.

        Attributes:
            id (int): Identificador único para cada cluster.
            pdb_reference_id (int): Clave foránea que referencia a la entidad 'PDBReference'. Se utiliza para
             identificar la estructura de proteína en la base de datos PDB a la que está asociado este cluster.
            chain_number (int): Número de cadenas de proteínas presentes en el cluster.
            cluster_id (String): Identificador del cluster, generalmente una cadena única que representa este grupo
             específico de cadenas de proteínas.
            is_representative (Boolean): Indica si el cluster es representativo de un conjunto más grande de cadenas
             similares. 'True' para sí, 'False' para no.
            sequence_length (int): Longitud promedio de las secuencias de las cadenas en el cluster.
            identity (Float): Valor que representa la identidad promedio de secuencia dentro del cluster,
             generalmente un porcentaje que indica cuán similares son las cadenas dentro del grupo.

        La relación con 'PDBReference' permite conectar cada cluster con su estructura correspondiente en la base de
         datos PDB, proporcionando un enlace directo a la información estructural detallada.
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
    Representa una cadena individual dentro de una estructura de proteína en la base de datos de estructuras de
     proteínas (PDB).

    Cada objeto `PDBChain` corresponde a una cadena polipeptídica específica en una estructura de proteína tal como se
     registra en PDB. Esta clase es crucial para el manejo detallado de estructuras de proteínas a nivel de cadena.

    Attributes:
        chains (String): Identificador de la cadena dentro de la estructura de proteína en PDB. Este campo es parte de la
         clave primaria compuesta.
        chain_number (String): Número o identificador adicional asociado a la cadena, formando parte de la clave
         primaria compuesta.
        sequence (String): Secuencia de aminoácidos de la cadena de proteína. Este campo es obligatorio y representa la
         secuencia lineal de aminoácidos de la cadena.
        pdb_reference_id (Integer): Clave foránea que referencia al identificador único de la estructura de proteína en
         la base de datos PDB. Este campo es parte de la clave primaria compuesta y establece una relación directa con
          la entidad `PDBReference`.
        pdb_reference (relationship): Relación con la entidad `PDBReference`, proporcionando un acceso directo a la
         información detallada de la estructura de proteína completa a la que pertenece esta cadena.

    La estructura de clave primaria compuesta de `chain`, `chain_number`, y `pdb_reference_id` asegura que cada
     instancia de `PDBChain` sea única y esté claramente vinculada a una estructura específica en PDB.
    """
    __tablename__ = 'pdb_chains'
    id = Column(Integer, primary_key=True)
    chains = Column(String)
    sequence = Column(String, nullable=False)
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'))

    pdb_reference = relationship("PDBReference", back_populates="pdb_chains")
    pdb_chain = relationship("PDBChain", back_populates="pdb_chains")



class PDBChain(Base):
    __tablename__ = 'pdb_chain'
    id = Column(Integer, primary_key=True)
    chain = Column(String)  # Ahora parte de la clave primaria
    # Foreign key to PDBChains
    pdb_chains_id = Column(Integer, ForeignKey('pdb_chains.id'))

    # Many-to-one relationship
    pdb_chains = relationship("PDBChains", back_populates="pdb_chain")



class UniProtPDBAlignment(Base):
    """
        Representa el alineamiento entre secuencias de proteínas de UniProt y la base de datos de estructuras de
         proteínas (PDB).

        Esta clase es crucial para correlacionar y comparar secuencias de proteínas entre las dos bases de datos más
         prominentes en bioinformática, UniProt y PDB. Permite a los investigadores y desarrolladores entender las
          similitudes y diferencias en la representación de secuencias de proteínas en diferentes bases de datos.

        Attributes:
            chain (String): Identificador de la cadena en el alineamiento. Este campo es parte de la clave primaria
             compuesta y puede referirse a una cadena específica en UniProt o PDB.
            pdb_reference_id (Integer): Clave foránea que referencia al identificador único de la estructura de proteína
             en PDB. Este campo es parte de la clave primaria compuesta y vincula el alineamiento a una estructura
              específica en PDB.
            uniprot_sequence_aligned (String): Secuencia de aminoácidos de UniProt que ha sido alineada con la secuencia
             de PDB.
            pdb_sequence_aligned (String): Secuencia de aminoácidos de PDB que ha sido alineada con la secuencia de
             UniProt.
            identity (Float): Valor numérico que representa el porcentaje de identidad entre las secuencias alineadas de
             UniProt y PDB.
            pdb_reference (relationship): Relación con la entidad `PDBReference`, proporcionando un enlace directo a la
             información estructural detallada de la proteína en PDB.

        La estructura de clave primaria compuesta de `chain` y `pdb_reference_id` asegura que cada instancia de
         `UniProtPDBAlignment` sea única y esté claramente vinculada a una estructura específica en PDB, facilitando el
          rastreo y análisis de alineaciones entre bases de datos.
    """
    __tablename__ = 'pdb_uniprot_chain_alignment'
    id = Column(Integer, primary_key=True)
    chain = Column(String)  # Ahora parte de la clave primaria
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'))  # Ahora parte de la clave primaria
    uniprot_sequence_aligned = Column(String)
    pdb_sequence_aligned = Column(String)
    identity = Column(Float)
    pdb_reference = relationship("PDBReference", back_populates="uniprot_pdb_alignments")


class CEAlignResults(Base):
    __tablename__ = 'ce_align_results'
    id = Column(Integer, primary_key=True)
    cluster_entry_id = Column(Integer, ForeignKey('clusters.id'))
    rms = Column(Float)


class GOTerm(Base):
    """
    Represents a Gene Ontology term associated with a protein.

    Attributes:
        id (int): Unique identifier for the GO term.
        go_id (str): Unique identifier in Gene Ontology.
        protein_entry_name (str): Entry name of the associated protein.
        protein (relationship): Relationship with the 'Protein' entity.
        category (str): Category of the GO term.
        description (str): Description of the GO term.
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

