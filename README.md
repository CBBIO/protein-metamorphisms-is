[![codecov](https://codecov.io/gh/CBBIO/protein-metamorphisms-is/graph/badge.svg?token=mtOqdG0xbU)](https://codecov.io/gh/CBBIO/protein-metamorphisms-is)
[![PyPI - Version](https://img.shields.io/pypi/v/protein-metamorphisms-is)](https://pypi.org/project/protein-metamorphisms-is/)
[![Documentation Status](https://readthedocs.org/projects/protein-metamorphisms-is/badge/?version=latest)](https://protein-metamorphisms-is.readthedocs.io/en/latest/?badge=latest)
![Linting Status](https://github.com/CBBIO/protein-metamorphisms-is/actions/workflows/test-lint.yml/badge.svg?branch=main)



# Exploration of “Metamorphism” and “Multi-functionality” in Proteins

💡 This study focuses on exploring phenomena of metamorphism and multifunctionality in proteins, fundamental aspects for understanding protein evolution and functionality across various biological contexts. We begin with a massive search for protein sequences that exhibit high percentages of identity, indicative of functional conservation across different species. Subsequently, we identify structures that, in addition to meeting this high identity criterion, exhibit significant differences in their spatial configuration, suggesting possible structural metamorphisms.

The main objective is to develop a comprehensive dataset that includes sets of varied structural conformations, providing a solid basis for comparative and evolutionary structural analysis. We implement an initial clustering strategy using the CD-HIT algorithm, which groups sequences based on their similarity. The resulting groups are then re-clustered using a structural embedding generation model based on the ID3 alphabet. This clustering focuses on finding within a homologous set of proteins those with greater structural differences, allowing a detailed analysis of their three-dimensional differences.

We have created an environment for analyzing structural distances that enables the analysis of inter-subcluster distances using leading structural alignment and distance algorithms such as Fatcat, US-Align, and CE-Align. High differences in these values are indicative of metamorphism. Concurrently, we employ similarity analysis techniques in Gene Ontology (GO) ontologies to discover proteins that, in addition to their structural conservation, exhibit multifunctionality. This includes an analysis of semantic distances per protein to identify disparate terms indicating multifunctionality.


> 📢 **FANTASIA Update in Progress**: 
> 
> An **updated version** of **FANTASIA** is currently being developed. This new version is a pipeline for **annotating GO (Gene Ontology) terms** in protein sequence files (FASTAs). Stay tuned for more details on the update.
> 
> This tool is being developed with **long-term support** in mind, and efforts are also underway to document its usage for **deployment on HPC** (High-Performance Computing) systems. Additionally, **FANTASIA** has been **re-engineered** to utilize updated dependencies and models under a formal development framework.
> 
> 🚀 **FANTASIA Support**: This framework has been **re-engineered** to use updated dependencies and models within a formal development framework. For more information, visit the [original project](https://github.com/MetazoaPhylogenomicsLab/FANTASIA).
> 
> 🔥 **Updated Documentation**: The most up-to-date and detailed documentation for this new version of **FANTASIA** is available [here](protein_metamorphisms_is/pipelines/fantasia/README.md).

## prerequisites

- Python 3.11.6
- RabbitMQ
- PostgreSQL with pgVector extension installed.


---

## Setup Instructions

### 1. Install Docker
Ensure Docker is installed on your system. If it’s not, you can download it from [here](https://docs.docker.com/get-docker/).

### 2. Set Up PostgreSQL with pgvector

Run the following command to start a PostgreSQL container with the pgvector extension:

```bash
docker run -d --name pgvectorsql \
    -e POSTGRES_USER=usuario \
    -e POSTGRES_PASSWORD=clave \
    -e POSTGRES_DB=BioData \
    -p 5432:5432 \
    pgvector/pgvector:pg16
```

Once the container is running, connect to the database and enable the `vector` extension:

```bash
docker exec -it pgvectorsql psql -U usuario -d BioData -c "CREATE EXTENSION IF NOT EXISTS vector;"
```

### 3. (Optional) Connect to the Database

You can use **pgAdmin 4**, a graphical interface for managing and interacting with PostgreSQL databases, or any other SQL client. Use the connection details defined during the `docker run` command (`POSTGRES_USER`, `POSTGRES_PASSWORD`, `POSTGRES_DB`, `port`) to access the database.

For more information or to download pgAdmin 4, visit the official website: [pgAdmin 4](https://www.pgadmin.org/).

After connecting to the database, you can:

- **Verify the database status**:
  Check that the database is running and accepting connections.
- **Inspect tables and schema**:
  Explore the database structure, including tables, columns, and relationships.
- **Run queries**:
  Perform SQL operations such as inserting, updating, or retrieving data.
- **Monitor activity**:
  Trace database activity by reviewing logs and active processes.

### 4. Set Up RabbitMQ

Start a RabbitMQ container using the command below:

```bash
docker run -d --name rabbitmq \
    -p 15672:15672 \
    -p 5672:5672 \
    rabbitmq:management
```

### 5. (Optional) Manage RabbitMQ

Once RabbitMQ is running, you can access its management interface to monitor and manage queues, exchanges, and messages.

Open your browser and go to: [RabbitMQ Management Interface](http://localhost:15672/#/queues).

From this interface, you can:

- **Inspect queues**:
  View the list of queues, their message rates, and other details.
- **Monitor message flow**:
  Track the rate of incoming and outgoing messages in real-time.

Ensure that the RabbitMQ container is running and accessible at `localhost:15672`. Use the default credentials (`guest`/`guest`) unless you have configured different ones.



--------
### Configuration Parameters and Constants

#### **Configuration Parameters (`config.yaml`)**

The main `config.yaml` file defines various settings used across the application, including system parameters, database configuration, and task-specific settings, pipelines such as FANTASIA contains a reduced template instance of this file on its own directory.

- **System Configuration**:
  - `max_workers`: Specifies the maximum number of worker threads that the system can use. This controls the concurrency level for processing tasks.
  - `binaries_path`: Defines the relative path where binary files required by the application are stored.

- **Database Configuration**:
  - `DB_USERNAME`: The username used to connect to the PostgreSQL database.
  - `DB_PASSWORD`: The password associated with the PostgreSQL user.
  - `DB_HOST`: The hostname of the PostgreSQL server (e.g., `localhost`).
  - `DB_PORT`: The port on which the PostgreSQL server is running (default: `5432`).
  - `DB_NAME`: The name of the PostgreSQL database where data is stored.

- **RabbitMQ Configuration**:
  - `rabbitmq_host`: The hostname of the RabbitMQ server.
  - `rabbitmq_user`: The username used to authenticate with RabbitMQ.
  - `rabbitmq_password`: The password associated with the RabbitMQ user.

- **Information System**:
  - `load_accesion_csv`: Path to the CSV file containing accession codes to be loaded into the system.
  - `load_accesion_column`: The column name in the CSV file that contains the accession IDs.
  - `tag`: A tag to be associated with the accession data for identification purposes.
  - `max_attempts`: Maximum number of attempts the system should make to process a given task before failing.
  - `search_criteria`: Criteria used to search for relevant data, typically in the format of a query string.
  - `limit`: The maximum number of results to return from a search query.

- **PDB Extraction**:
  - `resolution_threshold`: Maximum resolution (in Ångströms) for selecting PDB structures.
  - `server`: URL of the PDB server to download structure files.
  - `data_directory`: Path where PDB data files are stored.
  - `file_format`: Format of the PDB files (e.g., `mmCif`).
  - `allow_multiple_chain_models`: Boolean flag indicating whether to allow multiple chain models, typically set to `False` for NMR data.

- **Operations**:
  - `constants`: Path to the `constants.yaml` file containing various constant definitions used across the system.

- **Sequence Clustering**:
  - `fasta_path`: Path where the FASTA file containing sequences is saved.
  - `cdhit_out_path`: Path where CD-HIT output files are stored.
  - `sequence_identity_threshold`: Threshold for sequence identity used by CD-HIT to determine clusters.
  - `alignment_coverage`: Minimum alignment coverage required for sequences to be considered similar.
  - `memory_usage`: Maximum amount of memory (in MB) that CD-HIT is allowed to use.
  - `most_representative_search`: Boolean flag to enable or disable searching for the most representative sequence in each cluster.

- **Structural Alignment**:
  - `structural_alignment`: Contains settings for structural alignment tasks.
    - `types`: List of integers representing the types of alignments to perform.
    - `retry_timeout`: Time (in seconds) to wait before retrying a failed task.
    - `retry_count`: Maximum number of retry attempts for a failed task.
    - `batch_size`: Number of structures to process in a single batch.
    - `task_timeout`: Time (in seconds) to allow for the completion of a task before it times out.

- **Embedding**:
  - `embedding`: Contains settings for generating embeddings from sequences.
    - `types`: List of integers representing the types of embeddings to generate (e.g., ESM, Prost).
    - `batch_size`: Number of sequences to process in a single batch.

- **GO Metrics**:
  - `obo`: Path to the GO ontology file in OBO format.
  - `go_annotation_file`: Path to the file containing GO annotations.
  - `allowed_evidences`: List of evidence codes allowed when processing GO annotations.
  - `k`: Number of nearest neighbors to consider when predicting GO terms.

---

#### **Constants (`constants.yaml`)**

The `constants.yaml` file defines constant values used throughout the system, including types of structural alignments, complexity levels, embedding types, and prediction methods.

- **Structural Alignment Types**:
  - `name`: The name of the alignment method (e.g., "CE-align", "US-align").
  - `description`: A detailed explanation of the alignment method and its applications.
  - `task_name`: The internal name used by the system to reference the alignment task.

- **Structural Complexity Levels**:
  - `name`: The name representing the level of structural complexity (e.g., "Protein Chains", "Secondary Structures").
  - `description`: A description of what each complexity level represents in terms of protein structure.

- **Sequence Embedding Types**:
  - `name`: The name of the embedding model (e.g., "ESM", "Prost-T5").
  - `description`: A brief overview of the embedding model and its intended use.
  - `task_name`: The internal name used to reference the embedding task.
  - `model_name`: The specific model identifier used for generating embeddings.

- **Structure Embedding Types**:
  - `name`: The name of the structural embedding type (e.g., "3di").
  - `description`: A description of how the embedding type captures 3D structural information.
  - `task_name`: The internal name used for the embedding task.
  - `model_name`: The specific model identifier used for generating structural embeddings.

- **Prediction Methods**:
  - `name`: The name of the prediction method (e.g., "Cosine Similarity").
  - `description`: A description of how the prediction method works, particularly how it measures similarity between embeddings.

---


## **Get started:**

To execute the full process chain, simply run:

```bash
python main.py
```

This command will trigger the complete workflow, starting from the initial data preprocessing stages and continuing through to the final analysis and output generation.

## **Customizing the Workflow:**

You can customize the sequence of tasks executed by modifying `main.py` or adjusting the relevant parameters in the `config.yaml` file. This allows you to tailor the process flow to meet specific research needs or to experiment with different methodologies.
