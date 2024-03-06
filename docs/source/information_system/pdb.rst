PDB Module
==========

The PDB module is intricately designed for detailed interactions with the Protein Data Bank (PDB). It excels in downloading, processing, and managing complex 3D structural data of proteins, especially adept at handling PDB files containing multiple chains and models, using Biopython for efficient data manipulation.

Class Overview
--------------

The `PDBExtractor` class, extending `BioinfoExtractorBase`, is central to this module. It orchestrates the fetching and processing of 3D structural protein data from PDB, leveraging Biopython's capabilities for parsing and handling PDB files.

.. autoclass:: protein_metamorphisms_is.information_system.pdb.PDBExtractor
   :members:
   :undoc-members:
   :show-inheritance:

Key Features
------------

- **3D Structure Downloading**: Downloads protein 3D structures from PDB using Biopython's PDBList.
- **Data Processing with Biopython**: Utilizes Biopython's MMCIFParser and PDBIO for parsing PDB files and extracting chain information.
- **SQL/ORM Integration**: Seamlessly integrates with SQLAlchemy for storing and managing data in a relational database.
- **Concurrent Processing**: Employs multi-threading for efficient downloading and processing of multiple PDB structures.

Biopython Integration
---------------------

The module heavily relies on Biopython for various functionalities:

- **Chain Extraction**: Uses `ChainSelect <https://biopython.org/docs/latest/api/Bio.PDB.Dice.html?highlight=chainselect#Bio.PDB.Dice.ChainSelector>`_, a subclass of Biopython's Select, for extracting specific chains from PDB structures.
- **PDB File Handling**: Employs `PDBList <https://biopython.org/docs/latest/api/Bio.PDB.PDBList.html?highlight=pdb#module-Bio.PDB.PDBList>`_ for retrieving PDB files and MMCIFParser for parsing them.
- **Chain File Creation**: Utilizes `PDBIO <https://biopython.org/docs/latest/api/Bio.PDB.PDBIO.html?highlight=pdbio#module-Bio.PDB.PDBIO>`_ to write individual chain files for further analysis or storage.

SQL/ORM Entities
----------------

The module interacts with several SQL/ORM entities:

- :class:`~protein_metamorphisms_is.sql.model.PDBReference`: Manages references to Protein Data Bank entries.
- :class:`~protein_metamorphisms_is.sql.model.PDBChains`: Represents individual chains within a protein structure in the PDB.

These entities are crucial for organizing and storing the data fetched from the PDB in a relational database.

Configuration
-------------

The `PDBExtractor` class requires specific configuration settings:

.. code-block:: yaml

    # PDB Extraction Settings
    resolution_threshold: [Resolution threshold for PDB files]
    server: [PDB file server URL]
    pdb_path: [Local path to store PDB files]
    pdb_chains_path: [Local path to store PDB chains]
    file_format: "mmCif"  # The format of the PDB files. Currently, only the "mmCif" format is supported.


Adjust these settings based on your project requirements and available resources.

Usage
-----

To initiate the PDB data extraction process:

.. code-block:: python

   from protein_metamorphisms_is.information_system.pdb import PDBExtractor

   # Initialize the extractor with configuration
   pdb_extractor = PDBExtractor(conf)

   # Start the extraction process
   pdb_extractor.start()

This starts the extraction process, managing the download, processing, and storage of PDB data.
