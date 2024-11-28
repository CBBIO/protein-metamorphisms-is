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

