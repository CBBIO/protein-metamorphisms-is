from sqlalchemy import Column, Integer, String, Date, ForeignKey
from sqlalchemy.orm import relationship, declarative_base

Base = declarative_base()


class Proteina(Base):
    """
        Representa una proteína con sus propiedades y relaciones.

        Attributes:
            entry_name (str): Nombre de entrada único para la proteína, sirve
            como clave primaria.
            data_class (str): Clase de datos de la proteína.
            molecule_type (str): Tipo de molécula de la proteína.
            sequence_length (int): Longitud de la secuencia de la proteína.
            accessions (relationship): Relación con la entidad 'Accession'.
            created_date (Date): Fecha de creación del registro de la proteína.
            sequence_update_date (Date): Fecha de la última actualización
            de la secuencia.
            annotation_update_date (Date): Fecha de la última actualización de
             la anotación.
            description (str): Descripción de la proteína.
            gene_name (str): Nombre del gen asociado a la proteína.
            organism (str): Organismo del que proviene la proteína.
            organelle (str): Orgánulo asociado a la proteína.
            organism_classification (str): Clasificación del organismo.
            taxonomy_id (str): ID de taxonomía del organismo.
            host_organism (str): Organismo huésped de la proteína.
            host_taxonomy_id (str): ID de taxonomía del organismo huésped.
            comments (str): Comentarios adicionales sobre la proteína.
            pdb_references (relationship): Relación con la entidad
            'PDBReference'.
            go_terms (relationship): Relación con la entidad 'GOTerm'.
            keywords (str): Palabras clave asociadas a la proteína.
            protein_existence (int): Indicador de la existencia de la proteína.
            seqinfo (str): Información adicional sobre la secuencia.
        """
    __tablename__ = "proteinas"
    entry_name = Column(String, primary_key=True, unique=True, nullable=False)
    data_class = Column(String)
    molecule_type = Column(String)
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


class Accession(Base):
    """
        Representa un código de acceso para una proteína.

        Attributes:
            id (int): Identificador único para el acceso.
            accession_code (str): Código de acceso único para la proteína.
            proteina_entry_name (str): Nombre de entrada de la proteína
             asociada.
            proteina (relationship): Relación con la entidad 'Proteina'.
        """
    __tablename__ = "accessions"
    id = Column(Integer, primary_key=True)
    accession_code = Column(String, unique=True, nullable=False)
    proteina_entry_name = Column(String, ForeignKey("proteinas.entry_name"))
    proteina = relationship("Proteina", back_populates="accessions")


class PDBReference(Base):
    """
        Representa una referencia a la base de datos de estructuras de
        proteínas (PDB).

        Attributes:
            id (int): Identificador único para la referencia PDB.
            pdb_id (str): Identificador único en la base de datos PDB.
            proteina_entry_name (str): Nombre de entrada de la proteína
            asociada.
            proteina (relationship): Relación con la entidad 'Proteina'.
            method (str): Método utilizado para la determinación de la
            estructura.
            resolution (str): Resolución de la estructura.
            chains (str): Cadenas involucradas en la estructura.
        """
    __tablename__ = "pdb_references"
    id = Column(Integer, primary_key=True)
    pdb_id = Column(String, nullable=False)
    proteina_entry_name = Column(String, ForeignKey("proteinas.entry_name"))
    proteina = relationship("Proteina", back_populates="pdb_references")

    # Método utilizado para la determinación de la estructura
    method = Column(String)
    resolution = Column(String)  # Resolución de la estructura
    chains = Column(String)  # Cadenas involucradas


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
