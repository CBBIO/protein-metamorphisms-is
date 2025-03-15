Protein Information System (PIS)
========================================

Protein Information System (PIS) is an integrated biological information system focused on extracting, processing, and managing protein-related data. PIS consolidates data from UniProt, PDB, and GOA, enabling the efficient retrieval and organization of protein sequences, structures, and functional annotations.

The primary goal of PIS is to provide a robust framework for large-scale protein data extraction, facilitating downstream functional analysis and annotation transfer. The system is designed for high-performance computing (HPC) environments, ensuring scalability and efficiency.

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


Pipelines
===============
.. toctree::
   :maxdepth: 3
   :caption: Pipelines
   :hidden:


- `FANTASIA <https://github.com/CBBIO/FANTASIA>`_: Functional ANnoTAtion based on embedding space SImilArity.
    Using the information system, the FANTASIA pipeline has been set up and is now available in its corresponding repository.
    For full documentation you can visit `FANTASIA Documentation <https://fantasia.readthedocs.io/en/latest/>`_

- Massive search of metamorphisms and multifunctionality (WIP)


Extraction
===============
.. toctree::
   :maxdepth: 3
   :caption: Extraction
   :hidden:

   operation/extraction/accessions
   operation/extraction/uniprot

- `Accession Management <operation/extraction/accessions.html>`_: Manages the loading and processing of biological accession codes, ensuring data is correctly organized for subsequent analysis.
- `UniProt Extraction <operation/extraction/uniprot.html>`_: Handles the downloading and processing of detailed protein information from UniProt, enriching the database with accurate and up-to-date data.


Embedding
===============
.. toctree::
   :maxdepth: 3
   :caption: Embedding
   :hidden:

   operation/embedding/sequence_embedding

- `Sequence Embedding <operation/embedding/sequence_embedding.html>`_: Facilitates encoding of protein and nucleotide sequences into dense vectors.


Clustering
===============
.. toctree::
   :maxdepth: 3
   :caption: Clustering
   :hidden:

   operation/clustering/sequence_clustering
   operation/clustering/structural_subclustering

- `Sequence Clustering <operation/clustering/sequence_clustering.html>`_: Create different clusters using the CD-HIT algorithm based on sequence identity in amino acids.
- `Structural Subclustering <operation/clustering/structural_subclustering.html>`_: Within each cluster, create subclusters based on structural identity using structure embeddings.


Structural Alignment
=====================
.. toctree::
   :maxdepth: 3
   :caption: Structural Alignment
   :hidden:

   operation/structural_alignment/structural_alignment

- `Structural Alignment <operation/structural_alignment/structural_alignment.html>`_: Align the representative of each subcluster structurally.


Multifunctionality Analysis
=============================
.. toctree::
   :maxdepth: 3
   :caption: Multifunctionality Analysis
   :hidden:

   operation/functional/multifunctionality

- `Multifunctionality <operation/functional/multifunctionality.html>`_: Analyze semantic distances in Gene Ontology (GO) terms to detect multifunctionality.


FAQ
===============
.. toctree::
   :maxdepth: 2
   :caption:  Guides

   deployment/setup_instructions
