CDHit Module
============

The CDHit module in the Protein Data Handler project focuses on clustering protein sequences using the CD-HIT algorithm. It extends BioinfoOperatorBase for efficient database and configuration handling, catering to the needs of large-scale protein data processing.

Class Overview
--------------

The `CDHit` class is the centerpiece of this module, designed to handle the complexities of clustering protein sequences with CD-HIT. It manages loading sequences from the database, creating FASTA files, and executing the clustering algorithm.

.. autoclass:: protein_data_handler.operations.cdhit.CDHit
   :members:
   :undoc-members:
   :show-inheritance:

Key Features
------------

- **Sequence Clustering**: Efficiently clusters protein sequences using the CD-HIT algorithm.
- **Database Integration**: Seamlessly integrates with SQL databases for loading and storing protein data.
- **FASTA File Handling**: Capable of creating and managing FASTA files for the clustering process.
- **Multithreaded Operations**: Supports concurrent processing for improved performance.
- **Comprehensive Logging**: Implements robust logging for tracking the clustering process and troubleshooting.

Algorithm Integration
---------------------

`CDHit` utilizes the CD-HIT algorithm for its core functionality:

- **Efficient Clustering**: Employs the renowned CD-HIT algorithm for its efficiency in handling large sequence datasets.
- **Configurable Parameters**: Offers flexibility in setting CD-HIT parameters like sequence identity threshold and memory usage.

Installation of CD-HIT
----------------------

Before using the `CDHit` class, the CD-HIT algorithm must be installed on your system. CD-HIT can be easily installed on Debian-based systems using `apt-get`:

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install cd-hit

This ensures that the CD-HIT executable is available in your system's PATH, which is necessary for the `CDHit` class to function properly.

SQL/ORM Entities
----------------

The module interacts with several SQL/ORM entities for organizing and storing clustering results:

- :class:`~protein_data_handler.sql.model.PDBChains`: Manages Protein Data Bank chains information.
- :class:`~protein_data_handler.sql.model.Cluster`: Stores clustering results, including cluster identifiers and sequence information.

These entities facilitate the storage and retrieval of clustered protein sequence data in a relational database.

Configuration
-------------

`CDHit` requires specific configuration parameters for optimal operation. Here's a configuration template:

.. code-block:: yaml

    # CDHit Configuration
    fasta_path: [Path for FASTA file]
    cdhit_out_path: [Path for CD-HIT output file]
    sequence_identity_threshold: [Identity threshold for clustering]
    memory_usage: [Maximum memory usage for CD-HIT]
    max_workers: [Number of threads for CD-HIT]

Adjust these settings based on the specific requirements of your clustering tasks.

Usage
-----

To use `CDHit`, initialize the class with your configuration and start the clustering process:

.. code-block:: python

   from protein_data_handler.operations.cdhit import CDHit

   # Initialize with configuration
   cdhit_instance = CDHit(conf)

   # Start clustering
   cdhit_instance.start()

This initiates the sequence clustering process, handling all aspects from data loading to clustering and result storage.
