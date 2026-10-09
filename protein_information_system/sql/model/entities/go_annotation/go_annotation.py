from sqlalchemy import Column, ForeignKey, Integer, String, UniqueConstraint
from sqlalchemy.orm import relationship

from protein_information_system.sql.model.core.base import Base

""" 

Se ha modificado la clase para incluir el atributo publication en la tabla. 
Asimismo, se ha ajustado la definición de protein_id y go_id para que no formen 
parte de la clave primaria, evitando así la creación de una clave primaria 
compuesta (id + protein_id + go_id) que impediría registrar distintos códigos 
de evidencia para un mismo GO term.   

"""


class ProteinGOTermAnnotation(Base):
    __tablename__ = "protein_go_term_annotation"

    id = Column(
        Integer,
        primary_key=True,
        autoincrement=True,
    )

    protein_id = Column(
        String,
        ForeignKey("protein.id"),
        nullable=False,
        index=True,
    )

    go_id = Column(
        String,
        ForeignKey("go_terms.go_id"),
        nullable=False,
        index=True,
    )

    evidence_code = Column(
        String,
        nullable=False,
    )

    publication_id = Column(
        String,
        nullable=True,
        index=True,
    )

    __table_args__ = (
        UniqueConstraint(
            "protein_id",
            "go_id",
            "evidence_code",
            "publication_id",
            name="uq_pga_protein_go_evidence_publication",
        ),
    )

    protein = relationship(
        "Protein",
        back_populates="annotations",
    )

    go_term = relationship(
        "GOTerm",
        back_populates="annotations",
    )

    def __repr__(self):
        return (
            "<ProteinGOTermAnnotation("
            f"protein_id={self.protein_id}, "
            f"go_id={self.go_id}, "
            f"evidence_code={self.evidence_code}, "
            f"publication_id={self.publication_id}"
            ")>"
        )