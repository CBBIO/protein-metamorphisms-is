BioInformation System (BioInfo-IS)
==================================
The BioInformation System (BioInfo-IS) project, part of the computational suite in development, for exploring and
modeling different research areas related to **proteomics**, specially their structures and functionalities, it is
designed for the comprehensive management and processing of harmonized biological data.

BioInfo design (Plan)
==============================
- **BioInfo - This**: The core system that manages, processes, and organizes biological data, serving as the backbone of the BioResearch project.
- **BioInfo - Data Modeler**: A tool within the BioResearch suite that focuses on modeling and structuring biological data to enable advanced analysis and integration. It serves as the central hub where essential queries are executed, generating diverse datasets that are directly ready for processing by other tools.
- **BioInfo - Learning Models**: A component dedicated to machine learning and AI.
- **BioInfo - API**: Application Programming Interface that allows integration and access to system data and functionalities from other systems.
- **BioInfo - Portal**: Indicates an access platform or user interface where the project's tools and data are deployed for research exploitation.


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
   operation/extraction/pdb

- `Accession Management <operation/extraction/accessions.html>`_: Manages the loading and processing of biological accession codes, ensuring data is correctly organized for subsequent analysis.
- `UniProt Extraction <operation/extraction/uniprot.html>`_: Handles the downloading and processing of detailed protein information from UniProt, enriching the database with accurate and up-to-date data.
- `PDB Extraction <operation/extraction/pdb.html>`_:  Handles the extraction, processing, and storage of protein structures from the Protein Data Bank (PDB), managing the entire workflow from downloading files to storing processed data.


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


API
======
.. toctree::
   :maxdepth: 2
   :caption: API
   :hidden:

   tasks/base_api
   tasks/queue_api
   tasks/gpu_api

- `Base Api <tasks/base_api.html>`_: Esto es una prueba de la API
- `Queue Api <tasks/queue_api.html>`_: Otra API
- `GPU Api <tasks/queue_api.html>`_: La API de GPU


FAQ
===============
.. toctree::
   :maxdepth: 2
   :caption:  Guides

   deployment/setup_instructions
