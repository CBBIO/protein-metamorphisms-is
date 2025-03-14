[![codecov](https://codecov.io/gh/CBBIO/protein-metamorphisms-is/graph/badge.svg?token=mtOqdG0xbU)](https://codecov.io/gh/CBBIO/protein-metamorphisms-is)
[![PyPI - Version](https://img.shields.io/pypi/v/protein-metamorphisms-is)](https://pypi.org/project/protein-metamorphisms-is/)
[![Documentation Status](https://readthedocs.org/projects/protein-metamorphisms-is/badge/?version=latest)](https://protein-metamorphisms-is.readthedocs.io/en/latest/?badge=latest)
![Linting Status](https://github.com/CBBIO/protein-metamorphisms-is/actions/workflows/test-lint.yml/badge.svg?branch=main)

# **Protein Information System (PIS)**

**Protein Information System (PIS)** is an integrated biological information system focused on extracting, processing, and managing protein-related data. PIS consolidates data from **UniProt**, **PDB**, and **GOA**, enabling the efficient retrieval and organization of protein sequences, structures, and functional annotations.

The primary goal of PIS is to provide a robust framework for large-scale protein data extraction, facilitating downstream functional analysis and annotation transfer. The system is designed for **high-performance computing (HPC) environments**, ensuring scalability and efficiency.

## 📈 **Current State of the Project**

### **FANTASIA Redesign**
> 🔄 **FANTASIA has been completely redesigned and is now available at:**  
> [**FANTASIA Repository**](https://github.com/CBBIO/FANTASIA)  
> This new version is a pipeline for **annotating GO (Gene Ontology) terms** in protein sequence files (FASTAs). The redesign focuses on long-term support, updated dependencies, and improved integration with High-Performance Computing (HPC) environments.  

### **Stable Version of the Information System**
> 🛠️ **A stable version of the information system for working with UniProt and annotation transfer is available at:**  
> [**Zenodo Stable Release**](https://zenodo.org/records/14546346)  
> This version serves as a reference implementation and provides a consistent environment for annotation transfer tasks.

## **Prerequisites**

- Python 3.11.6
- RabbitMQ
- PostgreSQL with pgVector extension installed.

---

## **Setup Instructions**

### 1. Install Docker
Ensure Docker is installed on your system. If it’s not, you can download it from [here](https://docs.docker.com/get-docker/).

### 2. Starting Required Services

Ensure PostgreSQL and RabbitMQ services are running.

```bash
docker run -d --name pgvectorsql \
    --shm-size=64g \
    -e POSTGRES_USER=usuario \
    -e POSTGRES_PASSWORD=clave \
    -e POSTGRES_DB=BioData \
    -p 5432:5432 \
    pgvector/pgvector:pg16 \
    -c shared_buffers=16GB \
    -c effective_cache_size=32GB \
    -c work_mem=64MB
```

### 3. PostgreSQL Configuration

The configuration parameters provided above have been optimized for a machine with **128GB of RAM and 32 CPU cores**, allowing up to **20 concurrent workers**. These settings enhance PostgreSQL’s performance when handling large datasets and computationally intensive queries.

- `--shm-size=64g`: Allocates **64GB of shared memory** to the container, preventing PostgreSQL from running out of memory in high-performance environments.
- `-c shared_buffers=16GB`: Allocates **16GB of RAM** for PostgreSQL’s shared memory buffers. This should typically be **25-40%** of total system memory.
- `-c effective_cache_size=32GB`: Sets PostgreSQL’s estimated available memory for disk caching to **32GB**. This helps the query planner make better decisions.
- `-c work_mem=64MB`: Defines **64MB of memory per worker** for operations like sorting and hashing. This is crucial when handling parallel query execution.

### 4. (Optional) Connect to the Database

You can use **pgAdmin 4**, a graphical interface for managing and interacting with PostgreSQL databases, or any other SQL client.

### 5. Set Up RabbitMQ

Start a RabbitMQ container using the command below:

```bash
docker run -d --name rabbitmq \
    -p 15672:15672 \
    -p 5672:5672 \
    rabbitmq:management
```

### 6. (Optional) Manage RabbitMQ

Once RabbitMQ is running, you can access its management interface at [RabbitMQ Management Interface](http://localhost:15672/#/queues).

---

## **Get started:**

To execute the full extraction process, simply run:

```bash
python main.py
```

This command will trigger the complete workflow, starting from the initial data preprocessing stages and continuing through to the final data organization and storage.

## **Customizing the Workflow:**

You can customize the sequence of tasks executed by modifying `main.py` or adjusting the relevant parameters in the `config.yaml` file. This allows you to tailor the extraction process to meet specific research needs or to experiment with different data processing configurations.

