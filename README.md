[![codecov](https://codecov.io/gh/CBBIO/protein-metamorphisms-is/graph/badge.svg?token=mtOqdG0xbU)](https://codecov.io/gh/CBBIO/protein-metamorphisms-is)
[![PyPI - Version](https://img.shields.io/pypi/v/protein-metamorphisms-is)](https://pypi.org/project/protein-metamorphisms-is/)
[![Documentation Status](https://readthedocs.org/projects/protein-metamorphisms-is/badge/?version=latest)](https://protein-metamorphisms-is.readthedocs.io/en/latest/?badge=latest)
![Linting Status](https://github.com/CBBIO/protein-metamorphisms-is/actions/workflows/test-lint.yml/badge.svg?branch=main)

# Exploration of “Metamorphism” and “Multi-functionality” in Proteins

💡 This study focuses on exploring phenomena of metamorphism and multifunctionality in proteins, fundamental aspects for understanding protein evolution and functionality across various biological contexts. We begin with a massive search for protein sequences that exhibit high percentages of identity, indicative of functional conservation across different species. Subsequently, we identify structures that, in addition to meeting this high identity criterion, exhibit significant differences in their spatial configuration, suggesting possible structural metamorphisms.

The main objective is to develop a comprehensive dataset that includes sets of varied structural conformations, providing a solid basis for comparative and evolutionary structural analysis. We implement an initial clustering strategy using the CD-HIT algorithm, which groups sequences based on their similarity. The resulting groups are then re-clustered using a structural embedding generation model based on the ID3 alphabet. This clustering focuses on finding within a homologous set of proteins those with greater structural differences, allowing a detailed analysis of their three-dimensional differences.

We have created an environment for analyzing structural distances that enables the analysis of inter-subcluster distances using leading structural alignment and distance algorithms such as Fatcat, US-Align, and CE-Align. High differences in these values are indicative of metamorphism. Concurrently, we employ similarity analysis techniques in Gene Ontology (GO) ontologies to discover proteins that, in addition to their structural conservation, exhibit multifunctionality. This includes an analysis of semantic distances per protein to identify disparate terms indicating multifunctionality.

## 📈 **Current State of the Project**

### **FANTASIA Redesign**
> 🔄 **FANTASIA has been completely redesigned and is now available at:**  
> [**FANTASIA Repository**](https://github.com/CBBIO/FANTASIA)  
> This new version is a pipeline for **annotating GO (Gene Ontology) terms** in protein sequence files (FASTAs). The redesign focuses on long-term support, updated dependencies, and improved integration with High-Performance Computing (HPC) environments.  

### **Stable Version of the Information System**
> 🛠️ **A stable version of the information system for working with UniProt and annotation transfer is available at:**  
> [**Zenodo Stable Release**](https://zenodo.org/records/14546346)  
> This version serves as a reference implementation and provides a consistent environment for annotation transfer tasks.

## Prerequisites

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

You can use **pgAdmin 4**, a graphical interface for managing and interacting with PostgreSQL databases, or any other SQL client.

### 4. Set Up RabbitMQ

Start a RabbitMQ container using the command below:

```bash
docker run -d --name rabbitmq \
    -p 15672:15672 \
    -p 5672:5672 \
    rabbitmq:management
```

### 5. (Optional) Manage RabbitMQ

Once RabbitMQ is running, you can access its management interface at [RabbitMQ Management Interface](http://localhost:15672/#/queues).

---

## **Get started:**

To execute the full process chain, simply run:

```bash
python main.py
```

This command will trigger the complete workflow, starting from the initial data preprocessing stages and continuing through to the final analysis and output generation.

## **Customizing the Workflow:**

You can customize the sequence of tasks executed by modifying `main.py` or adjusting the relevant parameters in the `config.yaml` file. This allows you to tailor the process flow to meet specific research needs or to experiment with different methodologies.

