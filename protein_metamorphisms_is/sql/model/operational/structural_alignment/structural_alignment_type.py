from sqlalchemy import Column, Integer, String

from protein_metamorphisms_is.sql.model.core.base import Base


class StructuralAlignmentType(Base):
    """
    Provides a framework for aligning protein structures, crucial for understanding
    the functional and evolutionary relationships between proteins. This class
    enables the use of various alignment strategies, supporting a comprehensive
    approach to protein comparison.

    Structural alignment methods integrated within this framework include:

    - CE-align: Identifies optimal alignments based on the Combinatorial Extension
      method, focusing on similar backbone arrangements.
    - US-align: Utilizes an advanced algorithm for measuring structural similarity,
      offering insights into sequence identity and alignment scores.
    - FATCAT: Capable of accommodating protein flexibility during alignment,
      allowing for the detection of functionally important variations.

    By incorporating these methodologies, the class facilitates diverse approaches
    to protein comparison. This enables researchers to gain deeper insights into
    protein functionality and evolution, highlighting the significance of structural
    alignment in the field of bioinformatics.

    Attributes:
        id (Integer): Unique identifier for each alignment type.
        name (String): Name of the alignment type.
        description (String): Detailed description of the alignment method.
        task_name (String): Name of the specific task or process associated with this alignment type.
    """
    __tablename__ = 'structural_alignment_types'
    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)
    description = Column(String)
    task_name = Column(String)
