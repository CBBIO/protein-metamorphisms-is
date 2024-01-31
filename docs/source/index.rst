=============================================
Protein Data Handler: A Comprehensive Toolkit
=============================================

Welcome to the documentation of the Protein Data Handler, a versatile Python library designed to streamline the handling of protein data. This project offers a suite of tools for searching, downloading, and storing protein information from renowned databases like UniProt and the Protein Data Bank (PDB).

A distinctive feature of this library is its seamless integration with SQL databases, particularly PostgreSQL. It facilitates sophisticated querying and manipulation of data from UniProt and PDB at the protein chain level. Moreover, the Protein Data Handler supports the mass dumping of data from these databases into your own PostgreSQL database. This enables efficient management and analysis of large datasets, allowing users to perform complex SQL queries on protein data and execute in-depth analysis and operations on protein chains.

Whether you are a bioinformatician, a researcher in the field of proteomics, or just someone interested in protein data analysis, this library is crafted to cater to a wide range of needs. With an emphasis on ease of use and flexibility, the Protein Data Handler enables users to efficiently process and analyze protein-related data, with enhanced capabilities for robust data handling and customization.

Key Features
------------

- **UniProt Integration**: Seamlessly interact with the UniProt database for comprehensive protein data.
- **Protein Data Bank (PDB) Support**: Access and utilize data from the PDB, including 3D structures of proteins.
- **Database Management**: Efficiently manage and store protein data using robust database solutions.
- **User-Friendly Configuration**: Customize your experience with flexible YAML configuration options.

First Use Case: Exploring Metamorphisms
----------------------------------------

Our initial use case involves the exploration of metamorphisms in protein data. This process begins with the use of **CD-HIT**, an acclaimed bioinformatics tool for grouping similar sequence chains. CD-HIT aids in organizing biological sequences by reducing redundancy, thus enhancing processing speed and data manageability.

Once we have these groups formed through CD-HIT, we proceed with a detailed structural comparison using the **CE (Combinatorial Extension) alignment** method. This stage is vital as we select a representative sequence from each group and compare them against others to identify structural similarities and differences. The CE method proves highly effective in comparing protein structures, particularly in cases where sequences are similar, but their structures differ significantly. This approach is instrumental in our efforts to explore and understand metamorphisms in protein data.

Getting Started
---------------

Dive into the world of protein data handling by exploring the sections below:

.. toctree::
   :maxdepth: 2
   :caption: Information System

   uniprot
   pdb
   model

.. toctree::
   :maxdepth: 2
   :caption: Operations

   cdhit
   cealign

We hope this documentation serves as a valuable resource in your journey with Protein Data Handler. Happy data handling!
