database
========

El módulo `database` gestiona la conexión y la sesión con la base de datos utilizando SQLAlchemy.

Funciones
---------

.. autofunction:: protein_data_handler.helpers.database.database.create_session


Modelos ORM
-----------

Proteina
^^^^^^^^

.. autoclass:: protein_data_handler.models.uniprot.Proteina
   :members:
   :undoc-members:
   :show-inheritance:

Accession
^^^^^^^^^

.. autoclass:: protein_data_handler.models.uniprot.Accession
   :members:
   :undoc-members:
   :show-inheritance:

PDBReference
^^^^^^^^^^^^

.. autoclass:: protein_data_handler.models.uniprot.PDBReference
   :members:
   :undoc-members:
   :show-inheritance:

GOTerm
^^^^^^

.. autoclass:: protein_data_handler.models.uniprot.GOTerm
   :members:
   :undoc-members:
   :show-inheritance:
