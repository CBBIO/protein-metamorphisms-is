from sqlalchemy import Column, Integer, String, Float, ForeignKey
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base


class SequenceGoTermAnnotation(Base):
    __tablename__ = 'sequence_go_term_annotation'

    id = Column(Integer, primary_key=True, autoincrement=True)  # Clave primaria autoincremental
    go_id = Column(String, ForeignKey('go_terms.go_id'), nullable=False)  # Clave foránea con go_terms
    sequence_id = Column(Integer, ForeignKey('sequence.id'), nullable=False)  # Clave foránea con sequence
    distance = Column(Float, nullable=False)

    # Relaciones
    sequence = relationship("Sequence", back_populates="go_annotations")  # Relación con Sequence

    def __repr__(self):
        return f"<SequenceGOTermAnnotation(id={self.id}, go_id={self.go_id}, sequence_id={self.sequence_id}, distance={self.distance})>"
