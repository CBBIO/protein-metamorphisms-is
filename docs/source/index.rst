BioInformation System (BioInfo-IS)
==================================
The BioInformation System (BioInfo-IS) project, part of the BioResearch suite, is designed for the comprehensive management and processing of biological data, specifically in the field of computational biology. This project plays a critical role within BioResearch by managing tasks that handle data extraction, processing, and analysis, contributing to the broader goals of modeling and understanding biological processes.


BioInfo design
==============================
- **BioInfo - Information System**: The core system that manages, processes, and organizes biological data, serving as the backbone of the BioResearch project.
- **BioInfo - Data Modeler**: A tool within the BioResearch suite that focuses on modeling and structuring biological data to enable advanced analysis and integration. It serves as the central hub where essential queries are executed, generating diverse datasets that are directly ready for processing by other tools within the suite.
- **BioInfo - Learning Models**: A component dedicated to machine learning and AI.
- **BioInfo - API**: Application Programming Interface that allows integration and access to system data and functionalities from other systems.
- **BioInfo - Portal**: Indicates an access platform or user interface where the project's tools and data are centralized, offering a unified view.


Task Types
===============
.. toctree::
   :maxdepth: 2
   :caption: Task Types
   :hidden:

   tasks/base
   tasks/queue
   tasks/gpu

- `Base Tasks <tasks/base.html>`_: Handle Object Relational Model synchronization with the database and provide a session.
- `Queue Tasks <tasks/queue.html>`_: Handle queuing processes with RabbitMQ for task distribution and fault tolerance.
- `GPU Tasks <tasks/gpu.html>`_: Handle computationally intensive tasks using GPU.




Extraction Operations
===============
.. toctree::
   :maxdepth: 3
   :caption: Extraction
   :hidden:

   operation/extraction/accessions
   operation/extraction/uniprot
   operation/extraction/pdb

- `Accession Management <operation/extraction/accessions.html>`_: Manages the loading and processing of biological accession codes, ensuring data is correctly organized for subsequent analysis.
- `UniProt Extraction <operation/extraction/uniprot.html>`_: Handles the downloading and processing of detailed protein information from UniProt, enriching the database with accurate and up-to-date data.
- `PDB Extraction <operation/extraction/pdb.html>`_: Responsible for downloading and processing protein structures from the Protein Data Bank (PDB), allowing for the comparison and analysis of these structures within the system.

Clustering
----------

- `Sequence Clustering <operation/clustering/sequence_clustering.html>`_: Handles the clustering of protein sequences using the CD-HIT algorithm, reducing redundancy and organizing clustered data for subsequent analysis.

- `Sequence Based Clusters 3Di SubClustering <operation/clustering/sequence_embedding_subclustering.html>`_: Handles the subclustering of 3Di sequence embeddings within existing clusters, further refining the data for detailed structural analysis.

- `Optics Clustering <operation/optics_clustering.html>`_: Applies the OPTICS clustering algorithm to embeddings, identifying natural clusters within the data to facilitate further analysis.

Embedding
---------

- `Structure Embedding <operation/structure_embedding.html>`_: Generates 3Di sequences from protein structures, enabling their comparison and clustering based on structural configurations.

- `Sequence Embedding <operation/sequence_embedding.html>`_: Manages sequence embedding calculations using GPUs, optimizing computational resources for processing large volumes of data in an efficient manner.

Alignments
-------------

- `Multiple Structural Alignment (MSTA) <operation/multiple_structural_alignment.html>`_: Performs multiple structural alignments of proteins using FoldMason, comparing multiple structures simultaneously to identify similarities and differences.

- `RMSD Calculation <operation/rmsd_calculation.html>`_: Calculates the Root Mean Square Deviation (RMSD) of carbon atoms between aligned structures, providing a precise evaluation of structural similarities.

Funcional
---------

- `GO Prediction <operation/go_prediction.html>`_: Predicts Gene Ontology (GO) terms for protein sequences by identifying k-nearest neighbors in the embedding space and associating relevant GO terms.

Metrics
--------

- `GO Prediction Metrics <operation/go_prediction_metrics.html>`_: Evaluates the accuracy of GO term predictions by calculating semantic similarity metrics, such as Resnik similarity, between predicted and actual GO terms.
