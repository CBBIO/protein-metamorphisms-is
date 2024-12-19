from sqlalchemy import Column, Integer, ForeignKey, String, DateTime, func
from sqlalchemy.orm import relationship

from protein_metamorphisms_is.sql.model.core.base import Base


class Structure3Di(Base):
    __tablename__ = 'structure_3di'
    id = Column(Integer, primary_key=True)
    state_id = Column(Integer, ForeignKey('state.id'), nullable=False)
    embedding = Column(String)
    created_at = Column(DateTime, default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    state = relationship("State")
    subcluster_entries = relationship("SubclusterEntry", back_populates="structure_3di")
