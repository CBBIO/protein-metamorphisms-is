from sqlalchemy import Column, Integer, String, Float, ForeignKey
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base

class SequenceGoTermAnnotation(Base):
    __tablename__ = 'sequence_go_term_annotation'

    id = Column(Integer, primary_key=True)
    go_id = Column(String, ForeignKey('go_terms.go_id'), primary_key=True)
    embedding_id = Column(Integer, ForeignKey('sequence_embeddings.id'), unique=True)  # Relación 1:1 con SequenceEmbedding

    distance = Column(Float, nullable=False)
    embedding = relationship("SequenceEmbedding", back_populates="annotation")  # Relación con SequenceEmbedding

    def __repr__(self):
        return f"<SequenceGOTermAnnotation(protein_id={self.protein_id}, go_id={self.go_id}, evidence_code={self.evidence_code})>"
