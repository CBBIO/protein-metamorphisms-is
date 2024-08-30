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

## prerequisites

- python 3.10 or higher
- Access to a postgresql with pgVector extension installed.
```bash
docker run -d --name pgvectorsql \
    -e POSTGRES_USER=usuario \
    -e POSTGRES_PASSWORD=clave \
    -e POSTGRES_DB=BioData \
    -p 5432:5432 \
    pgvector/pgvector:pg16
```
- RabbitMQ
```bash

docker run -d --name rabbitmq \
    -p 15672:15672 \
    -p 5672:5672 \
    rabbitmq:management
```


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

## Get started 

the main file to start the system is `main.py`. you can run it as follows:

```sh
python main.py
```

## configuration files

- `config/config.yaml`: main configuration file containing system parameters, database configuration, and specific settings for each task.
- `config/constants.yaml`: defines constants for structural alignment types, levels of structural complexity, embedding types, and prediction methods.

