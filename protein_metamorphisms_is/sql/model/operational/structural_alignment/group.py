from sqlalchemy import Column, Integer, DateTime, func, ForeignKey
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base
from protein_metamorphisms_is.sql.model.operational.structural_alignment.structural_alignment_result import AlignmentResult


class AlignmentGroup(Base):
    __tablename__ = 'alignment_group'

    id = Column(Integer, primary_key=True, autoincrement=True)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    entries = relationship("AlignmentGroupEntry", back_populates="alignment_group", cascade="all, delete-orphan")
    result = relationship("AlignmentResult", back_populates="alignment_group", uselist=False)

    def __repr__(self):
        return f"<AlignmentGroup(id={self.id})>"



class AlignmentGroupEntry(Base):
    __tablename__ = 'alignment_group_entry'

    id = Column(Integer, primary_key=True, autoincrement=True)
    alignment_group_id = Column(Integer, ForeignKey('alignment_group.id'), nullable=False)
    cluster_entry_id = Column(Integer, ForeignKey('cluster_entry.id'), nullable=False)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    alignment_group = relationship("AlignmentGroup", back_populates="entries")
    cluster_entry = relationship("ClusterEntry")

    def __repr__(self):
        return f"<AlignmentGroupEntry(alignment_group_id={self.alignment_group_id}, cluster_entry_id={self.cluster_entry_id})>"
