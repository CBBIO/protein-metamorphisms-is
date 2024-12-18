from sqlalchemy import Column, DateTime, Boolean, Integer, String, func, Date, ForeignKey
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base


class Protein(Base):
    __tablename__ = "protein"

    id = Column(String, primary_key=True)
    sequence_id = Column(Integer, ForeignKey('sequence.id'))
    data_class = Column(String)
    molecule_type = Column(String)
    created_date = Column(Date)
    sequence_update_date = Column(Date)
    annotation_update_date = Column(Date)
    description = Column(String)
    gene_name = Column(String)
    organism = Column(String)
    organelle = Column(String)
    taxonomy_id = Column(String)
    comments = Column(String)
    protein_existence = Column(Integer)
    seqinfo = Column(String)
    disappeared = Column(Boolean)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    structure = relationship("Structure", back_populates="protein")
    sequence = relationship("Sequence", back_populates="protein", uselist=False)
    annotations = relationship("ProteinGOTermAnnotation", back_populates="protein")
    accessions = relationship("Accession", back_populates="protein")
    go_term_pairs = relationship("GOTermPair", back_populates="protein")  # Muchos go_term_pairs

    def __repr__(self):
        return f"<Protein(id={self.id}, gene_name={self.gene_name}, organism={self.organism})>"
