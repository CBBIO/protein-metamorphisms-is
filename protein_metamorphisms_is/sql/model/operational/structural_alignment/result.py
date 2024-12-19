from sqlalchemy import Column, Integer, Float, ForeignKey, DateTime, func
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base


class AlignmentResult(Base):
    __tablename__ = 'alignment_result'

    id = Column(Integer, primary_key=True, autoincrement=True)
    alignment_group_id = Column(Integer, ForeignKey('alignment_group.id'), nullable=False, unique=True)
    ce_rms = Column(Float, nullable=True)
    tm_rms = Column(Float, nullable=True)
    tm_seq_id = Column(Float, nullable=True)
    tm_score_chain_1 = Column(Float, nullable=True)
    tm_score_chain_2 = Column(Float, nullable=True)
    fc_rms = Column(Float, nullable=True)
    fc_identity = Column(Float, nullable=True)
    fc_similarity = Column(Float, nullable=True)
    fc_score = Column(Float, nullable=True)
    fc_align_len = Column(Float, nullable=True)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    alignment_group = relationship("AlignmentGroup", back_populates="result")

    def __repr__(self):
        return f"<AlignmentResult(alignment_group_id={self.alignment_group_id})>"
