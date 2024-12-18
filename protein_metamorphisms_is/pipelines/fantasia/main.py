"""
FANTASIA: Functional Annotation
========================================

FANTASIA (Functional ANnoTAtion based on embedding space SImilArity) is a pipeline for annotating GO terms in protein
sequence files using protein language models **ProtT5**, **ProstT5**, and **ESM2**.

FANTASIA takes as input a proteome file (either the longest isoform or the full set of isoforms for all genes),
removes identical sequences using **CD-HIT** and, optionally, excludes sequences longer than 5000 amino acids
(due to a length constraint in the model). It then computes embeddings for all sequences and stores them in an
**HDF5 file**. Finally, it queries a **vector database** to identify GO terms associated with similar proteins
based on the computed embeddings.

Pipeline Overview
-----------------
1. **Redundancy Filtering**:
   - Removes identical sequences using CD-HIT with a 95% similarity threshold.
   - Optionally excludes sequences longer than 5000 amino acids.

2. **Embedding Generation**:
   - Computes embeddings for protein sequences using **ProtT5**, **ProstT5**, and **ESM2**.
   - Stores embeddings in an HDF5 file, organized by sequence accession IDs and embedding types.

3. **GO Term Lookup**:
   - Embeddings are compared for similarity using a **vector database**.
   - Retrieves GO terms associated with the most similar proteins.
   - Results include GO terms, distances, and metadata.

4. **Results**:
   - Embeddings are saved in HDF5 files.
   - GO terms are saved in unique CSV files for each execution, named with a timestamp for reproducibility.

This pipeline results from joint efforts with equal contribution between **Ana Roja's lab**
(Andalusian Center for Developmental Biology, CSIC) and **Rosa Fernández's lab**
(Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, CSIC-UPF).
This collaboration highlights how synergistic efforts between labs with different expertise can achieve great outcomes.
We thank **LifeHUB-CSIC** for being the catalyst of this project and inspiring us to 'think big'.

Usage
-----
To run FANTASIA, ensure the `config.yaml` file is properly configured:

Example Execution:
------------------
.. code-block:: python

    python main.py --config /path/to/config.yaml

Key Configuration Parameters in `config.yaml`:
    - `fantasia_input_fasta`: Path to the input FASTA file.
    - `fantasia_output_h5`: Directory to store embeddings in HDF5 format.
    - `fantasia_output_csv`: Directory to store GO term results in CSV format.
    - `fantasia_prefix`: Prefix for naming output files.

Cite FANTASIA
-------------
Martínez-Redondo, G. I., Barrios, I., Vázquez-Valls, M., Rojas, A. M., & Fernández, R. (2024). Illuminating the
functional landscape of the dark proteome across the Animal Tree of Life.
https://doi.org/10.1101/2024.02.28.582465.

For our work about the performance of the different methods in model organisms, see:
Barrios-Núñez, I., Martínez-Redondo, G. I., Medina-Burgos, P., Cases, I., Fernández, R. & Rojas, A.M. (2024).
Decoding proteome functional information in model organisms using protein language models.
https://doi.org/10.1101/2024.02.14.580341.

Contact Information
-------------------
- Francisco Miguel Pérez Canales: fmpercan@upo.es
- Gemma I. Martínez-Redondo: gemma.martinez@ibe.upf-csic.es
- Ana M. Rojas: a.rojas.m@csic.es
- Rosa Fernández: rosa.fernandez@ibe.upf-csic.es
"""

from datetime import datetime
from protein_metamorphisms_is.pipelines.fantasia.embedder import SequenceEmbedder
from protein_metamorphisms_is.pipelines.fantasia.lookup import EmbeddingLookUp
from protein_metamorphisms_is.operation.extraction.accessions import AccessionManager
from protein_metamorphisms_is.operation.extraction.uniprot import UniProtExtractor
from protein_metamorphisms_is.operation.embedding.sequence_embedding import SequenceEmbeddingManager
from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config


def main(
        config_path="/home/bioxaxi/PycharmProjects/protein-metamorphisms-is/protein_metamorphisms_is/pipelines/fantasia/config.yaml"):
    conf = read_yaml_config(config_path)
    AccessionManager(conf).load_accessions_from_csv()
    UniProtExtractor(conf).start()
    SequenceEmbeddingManager(conf).start()
    current_date = datetime.now().strftime("%Y%m%d%H%M%S")

    # Inicializar procesos con current_date
    embedder = SequenceEmbedder(conf, current_date)
    embedder.start()

    lookup = EmbeddingLookUp(conf, current_date)
    lookup.start()


if __name__ == "__main__":
    main()
