from sqlalchemy import Column, Integer, DateTime, func, ForeignKey, String, Boolean
from sqlalchemy.orm import relationship
from protein_metamorphisms_is.sql.model.core.base import Base


class AlignmentGroup(Base):
    __tablename__ = 'alignment_group'

    id = Column(Integer, primary_key=True, autoincrement=True)
    comments = Column(String, default='No comments')
    is_metamorphic = Column(Boolean, default=False)
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
    subcluster_entry_id = Column(Integer, ForeignKey('subcluster_entry.id'), nullable=False, unique=False)  # Debe ser único para relación 1:1
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    # Relaciones
    alignment_group = relationship("AlignmentGroup", back_populates="entries")
    subcluster_entry = relationship("SubclusterEntry", back_populates="alignment_group_entry")  # Relación inversa uno a uno

    def __repr__(self):
        return f"<AlignmentGroupEntry(alignment_group_id={self.alignment_group_id}, subcluster_entry_id={self.subcluster_entry_id})>"
