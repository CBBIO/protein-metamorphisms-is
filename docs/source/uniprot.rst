UniProt Module
==============

The UniProt module in the Protein Data Handler project facilitates interactions with the UniProt database, focusing on extracting, processing, and storing protein data. It leverages SQL/ORM for efficient data management and integrates with a SQL database for persistence.

Class Overview
--------------

The `UniProtExtractor` class, an extension of `BioinfoExtractorBase`, is the core of this module. It manages the complexities of fetching and processing protein data from UniProt, including downloading records, parsing annotations, and integrating with a SQL database.

.. autoclass:: protein_data_handler.information_system.uniprot.UniProtExtractor
   :members:
   :undoc-members:
   :show-inheritance:

Key Features
------------

- **Data Extraction and Processing**: Fetches and processes protein data from UniProt.
- **Concurrent Downloads**: Uses multi-threading for efficient data retrieval.
- **SQL/ORM Integration**: Interacts with the database using SQLAlchemy for data storage and management.
- **Error Handling**: Implements error handling for reliable data extraction and processing.

Biopython Integration
---------------------

The `UniProtExtractor` utilizes Biopython for parsing and handling protein data:

- **ExPASy Access**: Uses `ExPASy <https://biopython.org/docs/latest/api/Bio.ExPASy.html>`_ to access and download protein information from the ExPASy database.
- **SwissProt Parsing**: Employs `SwissProt <https://biopython.org/docs/latest/api/Bio.SwissProt.html>`_ for parsing detailed protein information retrieved from UniProt.


SQL/ORM Entities
----------------

The module interacts with several SQL/ORM entities:

- :class:`~protein_data_handler.sql.model.Protein`: Unique identifiers for protein records in UniProt.
- :class:`~protein_data_handler.sql.model.Protein`: Stores detailed protein information.
- :class:`~protein_data_handler.sql.model.PDBReference`: Manages references to Protein Data Bank entries.
- :class:`~protein_data_handler.sql.model.UniprotChains`: Handles protein chain information linked to PDB entries.
- :class:`~protein_data_handler.sql.model.GOTerm`: Stores Gene Ontology terms associated with proteins.


These entities organize the data fetched from UniProt in a relational database.

Configuration
-------------

The `UniProtExtractor` class requires specific configuration settings for optimal operation. Below is a template for the configuration structure:

.. code-block:: yaml

    # System Configuration
    max_workers: [Number of concurrent workers]

    # Database Configuration
    DB_USERNAME: [Database username]
    DB_PASSWORD: [Database password]
    DB_HOST: [Database host address]
    DB_PORT: [Database port]
    DB_NAME: [Database name]

    # UniProt Extraction Settings
    search_criteria: [UniProt search criteria]
    limit: [Maximum number of records to fetch per page (does not affect)]

Adjust these settings based on your project requirements and available resources.

Usage
-----

Initialize `UniProtExtractor` with the configuration and start the extraction process:

.. code-block:: python

   from protein_data_handler.information_system.uniprot import UniProtExtractor

   # Initialize with configuration
   extractor = UniProtExtractor(conf)

   # Start extraction
   extractor.start()

This initiates the extraction process, managing the download, processing, and storage of UniProt data.

