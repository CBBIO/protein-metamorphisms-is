Fantasia Parameters
===================

The configuration file is divided into sections to control various aspects of the pipeline. Below is a summary of the key parameters:
Default system configuration should be OK for any Unix system, for the current state of GOA.

1. System Settings
------------------

- `max_workers`: Number of workers for parallel processing. Default: `10`.

2. Database Settings
--------------------

- `DB_USERNAME`: PostgreSQL username. Default: `usuario`.
- `DB_PASSWORD`: PostgreSQL password. Default: `clave`.
- `DB_HOST`: Host address for the database. Default: `localhost`.
- `DB_PORT`: Port for database connection. Default: `5432`.
- `DB_NAME`: Name of the database. Default: `BioData`.

3. RabbitMQ Settings
---------------------

- `rabbitmq_host`: Host for RabbitMQ. Default: `localhost`.
- `rabbitmq_user`: RabbitMQ username. Default: `guest`.
- `rabbitmq_password`: RabbitMQ password. Default: `guest`.

4. Execution Settings
----------------------

- `limit_execution`: Maximum number of items to process (for debugging). Default: `100`.

5. LookUp table from CSV
-------------------------

- `load_accesion_csv`: Path to the CSV file containing accessions. Default: `~/fantasia/input/genes.csv`.
- `load_accesion_column`: Column name in the CSV file. Default: `id`.
- `tag`: Identifier tag for accession loading. Default: `'GOA2022'`.

5b. LookUp table from UniProt (By default)
------------------------------------------------------------------

- `search_criteria`: Query string for UniProt search. Default: `'(reviewed:true)'`.
- `limit`: Number of results to fetch from UniProt. Default: `200`.
- `tag`: Identifier tag for accession loading. Default: `'GOA2022'`. (It's the same defined in 5.)

6. Evidence Codes
-----------------------------------------------------------------------------------
- `allowed_evidences`: This parameter specifies the list of UniProt evidence codes used to filter annotations. The default selection includes codes that indicate experimental evidence:

  - `EXP`: Inferred from Experiment.
  - `IDA`: Inferred from Direct Assay.
  - `IPI`: Inferred from Physical Interaction.
  - `IMP`: Inferred from Mutant Phenotype.
  - `IGI`: Inferred from Genetic Interaction.
  - `IEP`: Inferred from Expression Pattern.
  - `TAS`: Traceable Author Statement.
  - `IC`: Inferred by Curator.

  These default codes represent high-confidence annotations derived from experimental and curated evidence.

  However, you can include any evidence codes from the full set defined by the Gene Ontology Consortium. For more details on all available evidence codes, visit the `Guide to GO Evidence Codes <https://geneontology.org/docs/guide-go-evidence-codes/>`_.

7. FANTASIA Settings
-------------------------------------------

- `fantasia_input_fasta`: Path to the input FASTA file. Default: `~/fantasia/input/input.fasta`.
- `fantasia_output_h5`: Directory to save embeddings in HDF5 format. Default: `~/fantasia/embeddings/`.
- `fantasia_output_csv`: Directory to save GO term results in CSV format. Default: `~/fantasia/results/`.
- `fantasia_prefix`: Prefix for output files. Default: `finger_zinc`.
- `max_distance`: Maximum similarity distance for GO term matches. Default: `1`.
- `length_filter`: Maximum allowed sequence length. Default: `5000`.
- `redundancy_filter`: Threshold for redundancy removal. Default: `0.95`.
- `redundancy_file`: Path to save redundancy-filtered FASTA. Default: `~/fantasia/redundancy/output.fasta`.
- `topgo`: Enable or disable TopGO analysis. Default: `True`.

8. Embedding Settings
--------------------------------------------

- `embedding.types`: List of embedding types to compute (Comment with # at the beginning to de/activate):

  - `1`: ESM2.
  - `2`: ProstT5.
  - `3`: ProtT5.
- `batch_size`: Enqueue by this batch size. Default: `2000`
- `embedding.batch_size`: batch size to feed models. Default: `40`.

9. System Constants
--------------------

The primary constant value for the FANTASIA pipeline includes pre-configured models for sequence embedding. You can easily expand this list by adding your own models to the configuration file and integrating them into the pipeline with minimal effort.

- `constants`: Path to the constants configuration file. Default: `config/constants.yaml`.