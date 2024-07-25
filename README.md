[![codecov](https://codecov.io/gh/CBBIO/protein-metamorphisms-is/graph/badge.svg?token=mtOqdG0xbU)](https://codecov.io/gh/CBBIO/protein-metamorphisms-is)
[![PyPI - Version](https://img.shields.io/pypi/v/protein-metamorphisms-is)](https://pypi.org/project/protein-metamorphisms-is/)
[![Documentation Status](https://readthedocs.org/projects/protein-metamorphisms-is/badge/?version=latest)](https://protein-metamorphisms-is.readthedocs.io/es/latest/?badge=latest)
![Linting Status](https://github.com/CBBIO/protein-metamorphisms-is/actions/workflows/test-lint.yml/badge.svg?branch=main)



# Exploration of “Metamorphism” and “Multi-functionality” in Proteins

💡 This study focuses on exploring phenomena of metamorphism and multifunctionality in proteins, fundamental aspects for understanding protein evolution and functionality across various biological contexts. We begin with a massive search for protein sequences that exhibit high percentages of identity, indicative of functional conservation across different species. Subsequently, we identify structures that, in addition to meeting this high identity criterion, exhibit significant differences in their spatial configuration, suggesting possible structural metamorphisms.

The main objective is to develop a comprehensive dataset that includes sets of varied structural conformations, providing a solid basis for comparative and evolutionary structural analysis. We implement an initial clustering strategy using the CD-HIT algorithm, which groups sequences based on their similarity. The resulting groups are then re-clustered using a structural embedding generation model based on the ID3 alphabet. This clustering focuses on finding within a homologous set of proteins those with greater structural differences, allowing a detailed analysis of their three-dimensional differences.

We have created an environment for analyzing structural distances that enables the analysis of inter-subcluster distances using leading structural alignment and distance algorithms such as Fatcat, US-Align, and CE-Align. High differences in these values are indicative of metamorphism. Concurrently, we employ similarity analysis techniques in Gene Ontology (GO) ontologies to discover proteins that, in addition to their structural conservation, exhibit multifunctionality. This includes an analysis of semantic distances per protein to identify disparate terms indicating multifunctionality.

Additionally, we utilize sequence embeddings based on ProstT5, ProtT5, ESM models to perform automatic transfer of GO terms. This transfer is not only carried out through a single type of embeddings but also through the possible concatenation of these, thus enriching the precision of our predictions.

Furthermore, we are developing predictors that may indicate multifunctionality or metamorphism, through the filtering of the datasets resulting from these operations. This integrated approach not only expands our understanding of the structural and functional plasticity of proteins but also significantly contributes to bioinformatics and structural biology, providing insights into the adaptability and evolution of proteins over time.

## project structure

the project is structured as follows:

```markdown
protein-metamorphisms-is/
│
├── __init__.py
│
├── base/
│   ├── base.py
│   ├── gpu.py
│   └── queue.py
│
├── config/
│   ├── config.yaml
│   └── constants.yaml
│
├── helpers/
│   ├── __init__.py
│   ├── config/
│   │   ├── __init__.py
│   │   └── yaml.py
│   ├── logger/
│   │   ├── __init__.py
│   │   └── logger.py
│   └── parser/
│       ├── __init__.py
│       └── parser.py
│
├── information_system/
│   ├── __init__.py
│   ├── accessions.py
│   ├── pdb.py
│   └── uniprot.py
│
├── main.py
│
├── operations/
│   ├── __init__.py
│   ├── base/
│   │   ├── __init__.py
│   │   └── operator.py
│   ├── cdhit.py
│   ├── cuda.py
│   ├── embedding_tasks/
│   │   ├── __init__.py
│   │   └── esm.py
│   │   └── prost_t5.py
│   ├── go_metrics.py
│   ├── go_prediction.py
│   ├── optics.py
│   ├── protein_go_prediction_metrics.py
│   ├── seq_embeddings.py
│   ├── structural_alignment.py
│   ├── structural_alignment_tasks/
│   │   ├── __init__.py
│   │   ├── combinatorial_extension.py
│   │   ├── fatcat.py
│   │   └── universal.py
│   └── structure_embeddings.py
│
├── sql/
│   ├── __init__.py
│   ├── base/
│   │   ├── database_manager.py
│   ├── constants.py
│   └── model.py
```

## prerequisites

- python 3.10 or higher
- necessary python libraries (see `requirements.txt`)
- access to a postgresql database
- rabbitmq for queue management

## installation

clone the repository:

```sh
git clone https://github.com/CBBIO/protein-metamorphisms-is.git
cd protein-metamorphisms-is
```

install the dependencies:

```sh
poetry install
```

configure the database and rabbitmq by editing the `config/config.yaml` file:

```yaml

max_workers: 10
binaries_path: '../binaries'

db_username: user
db_password: password
db_host: localhost
db_port: 5432
db_name: biodata

rabbitmq_host: localhost
rabbitmq_user: guest
rabbitmq_password: guest
...
```

## usage

the main file to start the system is `main.py`. you can run it as follows:

```sh
python main.py
```

this file initializes various components based on the configuration and starts different processes for data extraction and processing.

## configuration structure

- `config/config.yaml`: main configuration file containing system parameters, database configuration, and specific settings for each task.
- `config/constants.yaml`: defines constants for structural alignment types, levels of structural complexity, embedding types, and prediction methods.

## main components

### base

- `base.py`: abstract base class for initializing tasks, including configuration loading, session initialization, and abstract methods for starting, processing, and storing data.
- `gpu.py`: class for initializing gpu-based tasks using rabbitmq for queuing and multiprocessing for parallel processing.
- `queue.py`: base class for queue-based tasks using rabbitmq, including setup, worker initialization, and message processing.

### helpers

- `logger/logger.py`: logger configuration for recording events and debug messages.

### information_system

- `accessions.py`: management of protein accessions.
- `pdb.py`: extraction of pdb data.
- `uniprot.py`: extraction of uniprot data.

### operations

- `cdhit.py`: management of sequence clustering using cd-hit.
- `cuda.py`: tasks related to cuda.
- `embedding_tasks/esm.py`: embedding tasks using esm.
- `embedding_tasks/prost_t5.py`: embedding tasks using prot-t5.
- `go_metrics.py`: calculation of go metrics.
- `go_prediction.py`: prediction of go terms.
- `seq_embeddings.py`: management of sequence embeddings.
- `structural_alignment.py`: management of structural alignments.
