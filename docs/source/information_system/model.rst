ORM Model
==========

This section provides detailed documentation of the Object-Relational Mapping (ORM) models used in the Protein Data Handler. Each model represents a table in the database and is crucial for managing and storing protein data efficiently.

ER Diagram
----------
.. image:: ../../images/model.png
   :alt: Model Diagram
   :align: center

Protein
-------

.. autoclass:: protein_metamorphisms_is.sql.model.Protein
   :members:
   :show-inheritance:

Accession
---------

.. autoclass:: protein_metamorphisms_is.sql.model.Accession
   :members:
   :show-inheritance:

PDBReference
------------

.. autoclass:: protein_metamorphisms_is.sql.model.PDBReference
   :members:
   :show-inheritance:

PDBChains
---------

.. autoclass:: protein_metamorphisms_is.sql.model.PDBChains
   :members:
   :show-inheritance:

Cluster
-------

.. autoclass:: protein_metamorphisms_is.sql.model.Cluster
   :members:
   :show-inheritance:

StructuralComplexityLevel
-------------------------

.. autoclass:: protein_metamorphisms_is.sql.model.StructuralComplexityLevel
   :members:
   :show-inheritance:

StructuralAlignmentType
-----------------------

.. autoclass:: protein_metamorphisms_is.sql.model.StructuralAlignmentType
   :members:
   :show-inheritance:

StructuralAlignmentQueue
------------------------

.. autoclass:: protein_metamorphisms_is.sql.model.StructuralAlignmentQueue
   :members:
   :show-inheritance:

StructuralAlignmentResults
--------------------------

.. autoclass:: protein_metamorphisms_is.sql.model.StructuralAlignmentResults
   :members:
   :show-inheritance:

GOTerm
------

.. autoclass:: protein_metamorphisms_is.sql.model.GOTerm
   :members:
   :show-inheritance:
