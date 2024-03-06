===========================================================
Information System for Protein Metamorphisms discovery.
===========================================================

Welcome to the documentation of the Protein Data Handler for Metamorphisms discovery, a versatile Python library designed to streamline the handling of protein data. This project offers a suite of tools for searching, downloading, and storing protein information from renowned databases like UniProt and the Protein Data Bank (PDB).

A distinctive feature of this library is its seamless integration with SQL databases, particularly PostgreSQL. It facilitates sophisticated querying and manipulation of data from UniProt and PDB at the protein chain level. Moreover, the Protein Data Handler supports the mass dumping of data from these databases into your own PostgreSQL database. This enables efficient management and analysis of large datasets, allowing users to perform complex SQL queries on protein data and execute in-depth analysis and operations on protein chains.

Our use case involves the exploration of metamorphisms in protein data. This process begins with the use of **CD-HIT**, an acclaimed bioinformatics tool for grouping similar sequence chains. CD-HIT aids in organizing biological sequences by reducing redundancy, thus enhancing processing speed and data manageability.

Once we have these groups formed through CD-HIT, we proceed with a detailed structural alignment comparison using multiple standard methods. This stage is vital as we select a representative sequence from each group and compare them against others to identify structural similarities and differences. Those methods proves highly effective in comparing protein structures, particularly in cases where sequences are similar, but their structures differ significantly. This approach is instrumental in our efforts to explore and understand metamorphisms in protein data.

Key Features
------------

- **UniProt Integration**: Seamlessly interact with the UniProt database for comprehensive protein data.
- **Protein Data Bank (PDB) Support**: Access and utilize data from the PDB, including 3D structures of proteins.
- **Database Management**: Efficiently manage and store protein data using robust database solutions. Enhance your database interactions with the integration of ORM (Object-Relational Mapping) for streamlined data operations and a well-structured SQL model for optimized data organization.
- **User-Friendly Configuration**: Customize your experience with flexible YAML configuration options.
- **Robust Computational Result Acquisition**: Ensure reliability in the processing of computational results through a system of queuing and concurrent processing.


Getting Started
---------------

Dive into the world of protein metamorphisms discovery by exploring the sections below:

.. toctree::
   :maxdepth: 3
   :caption: Information System

   information_system/uniprot
   information_system/pdb
   information_system/model

.. toctree::
   :maxdepth: 3
   :caption: Operations

   operation/cdhit
   operation/structural_alignment

We hope this documentation serves as a valuable resource in your journey with Protein Data Handler. Happy data handling!
