from sqlalchemy import (Column, Integer, String, Date, ForeignKey, DateTime,
                        func, Float, Boolean)
from sqlalchemy.orm import relationship, declarative_base

Base = declarative_base()


class Proteina(Base):
    """
    Representa una proteína con sus propiedades y relaciones en una base de datos.

    Attributes:
        entry_name (str): Nombre único de entrada para la proteína, utilizado como clave primaria.
        data_class (str): Clasificación de los datos de la proteína.
        molecule_type (str): Tipo de molécula de la proteína.
        sequence_length (int): Longitud de la secuencia de aminoácidos de la proteína.
        accessions (relationship): Relación con la clase 'Accession', representando los códigos de acceso asociados a la
         proteína.
        created_date (Date): Fecha en la que se creó el registro de la proteína.
        sequence_update_date (Date): Fecha de la última actualización de la secuencia de la proteína.
        annotation_update_date (Date): Fecha de la última actualización de la anotación de la proteína.
        description (str): Descripción general de la proteína.
        gene_name (str): Nombre del gen asociado a esta proteína.
        organism (str): Organismo al que pertenece la proteína.
        organelle (str): Orgánulo específico donde se encuentra la proteína, si aplica.
        organism_classification (str): Clasificación taxonómica del organismo.
        taxonomy_id (str): Identificador de taxonomía para el organismo.
        host_organism (str): Organismo huésped de la proteína, si aplica.
        host_taxonomy_id (str): Identificador de taxonomía para el organismo huésped.
        comments (str): Comentarios adicionales sobre la proteína.
        pdb_references (relationship): Relación con la clase 'PDBReference', que contiene referencias a la base de
         datos de estructuras de proteínas.
        go_terms (relationship): Relación con la clase 'GOTerm', representando términos de Gene Ontology asociados a la
         proteína.
        keywords (str): Palabras clave relacionadas con la proteína.
        protein_existence (int): Indicador numérico de la existencia de la proteína.
        seqinfo (str): Información adicional sobre la secuencia de la proteína.
        disappeared (Boolean): Indica si la proteína ya no está presente o relevante.
        created_at (DateTime): Fecha y hora en la que se creó el registro.
        updated_at (DateTime): Fecha y hora de la última actualización del registro.
    """
    __tablename__ = "proteinas"
    entry_name = Column(String, primary_key=True, unique=True, nullable=False)
    data_class = Column(String)
    molecule_type = Column(String)
    sequence = Column(String)
    sequence_length = Column(Integer)
    accessions = relationship(
        "Accession", back_populates="proteina"
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

    pdb_references = relationship("PDBReference", back_populates="proteina")
    go_terms = relationship("GOTerm", back_populates="proteina")

    keywords = Column(String)  # Similar a comments
    protein_existence = Column(Integer)
    seqinfo = Column(String)
    disappeared = Column(Boolean)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())


class Accession(Base):
    """
    Representa un código de acceso único para una proteína en una base de datos.

    Attributes:
        id (int): Identificador único para el acceso.
        accession_code (str): Código de acceso único para la proteína.
        proteina_entry_name (str): Nombre de entrada de la proteína asociada.
        proteina (relationship): Relación con la clase 'Proteina'.
        disappeared (Boolean): Indica si el código de acceso ya no está en uso.
        created_at (DateTime): Fecha y hora de creación del registro.
        updated_at (DateTime): Fecha y hora de la última actualización del registro.
    """
    __tablename__ = "accessions"
    id = Column(Integer, primary_key=True)
    accession_code = Column(String, unique=True, nullable=False)
    proteina_entry_name = Column(String, ForeignKey("proteinas.entry_name"))
    proteina = relationship("Proteina", back_populates="accessions")
    disappeared = Column(Boolean)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())


class PDBReference(Base):
    """
    Representa una referencia a la base de datos de estructuras de proteínas (PDB).

    Attributes:
        id (int): Identificador único para la referencia PDB.
        pdb_id (str): Identificador único en la base de datos PDB.
        proteina_entry_name (str): Nombre de entrada de la proteína asociada.
        proteina (relationship): Relación con la clase 'Proteina'.
        method (str): Método utilizado para la determinación de la estructura de la proteína.
        resolution (Float): Resolución de la estructura de la proteína en la base de datos PDB.
        uniprot_chains (relationship): Relación con la clase 'UniprotChain', representando las cadenas de Uniprot en la
         estructura de la proteína.
        pdb_chains (relationship): Relación con la clase 'PDBChain', representando las cadenas de la estructura en la
         base de datos PDB.
        uniprot_pdb_alignments (relationship): Relación con la clase 'UniProtPDBAlignment', representando alineaciones
         entre secuencias de Uniprot y PDB.
        created_at (DateTime): Fecha y hora de creación del registro.
        updated_at (DateTime): Fecha y hora de la última actualización del registro.
    """
    __tablename__ = "pdb_references"
    id = Column(Integer, primary_key=True)
    pdb_id = Column(String, nullable=False)
    proteina_entry_name = Column(String, ForeignKey("proteinas.entry_name"))
    proteina = relationship("Proteina", back_populates="pdb_references")

    # Método utilizado para la determinación de la estructura
    method = Column(String)
    resolution = Column(Float)  # Resolución de la estructura
    uniprot_chains = relationship("UniprotChain", back_populates="pdb_reference")
    pdb_chains = relationship("PDBChain", back_populates="pdb_reference")
    uniprot_pdb_alignments = relationship("UniProtPDBAlignment", back_populates="pdb_reference")
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())


class UniprotChain(Base):
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
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'),
                              primary_key=True)  # Ahora parte de la clave primaria
    chain = Column(String, primary_key=True)  # Ahora parte de la clave primaria
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
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'))
    chain_number = Column(Integer)
    cluster_id = Column(String)
    is_representative = Column(Boolean)  # 0: False, 1: True
    sequence_length = Column(Integer)
    identity = Column(Float)


class PDBChain(Base):
    """
    Representa una cadena individual dentro de una estructura de proteína en la base de datos de estructuras de
     proteínas (PDB).

    Cada objeto `PDBChain` corresponde a una cadena polipeptídica específica en una estructura de proteína tal como se
     registra en PDB. Esta clase es crucial para el manejo detallado de estructuras de proteínas a nivel de cadena.

    Attributes:
        chain (String): Identificador de la cadena dentro de la estructura de proteína en PDB. Este campo es parte de la
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

    chain = Column(String)  # Parte de la clave primaria
    chain_number = Column(String, primary_key=True)
    sequence = Column(String, nullable=False)
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'),
                              primary_key=True)  # Parte de la clave primaria

    pdb_reference = relationship("PDBReference", back_populates="pdb_chains")


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

    chain = Column(String, primary_key=True)  # Ahora parte de la clave primaria
    pdb_reference_id = Column(Integer, ForeignKey('pdb_references.id'),
                              primary_key=True)  # Ahora parte de la clave primaria
    uniprot_sequence_aligned = Column(String)
    pdb_sequence_aligned = Column(String)
    identity = Column(Float)
    pdb_reference = relationship("PDBReference", back_populates="uniprot_pdb_alignments")


class GOTerm(Base):
    """
        Representa un término del Gene Ontology asociado a una proteína.

        Attributes:
            id (int): Identificador único para el término GO.
            go_id (str): Identificador único en Gene Ontology.
            proteina_entry_name (str): Nombre de entrada de la
            proteína asociada.
            proteina (relationship): Relación con la entidad 'Proteina'.
            category (str): Categoría del término GO.
            description (str): Descripción del término GO.
        """
    __tablename__ = "go_terms"
    id = Column(Integer, primary_key=True)
    go_id = Column(String, nullable=False)
    proteina_entry_name = Column(String, ForeignKey("proteinas.entry_name"))
    proteina = relationship("Proteina", back_populates="go_terms")
    category = Column(String)
    description = Column(String)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())
