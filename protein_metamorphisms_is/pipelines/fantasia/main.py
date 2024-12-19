"""
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
   - Removes identical sequences using CD-HIT with a 95% similarity threshold. (Optional)
   - Optionally excludes sequences longer than 5000 amino acids. (Optional)

2. **Embedding Generation**:
   - Computes embeddings for protein sequences using **ProtT5**, **ProstT5**, and **ESM2**.
   - Stores embeddings in an HDF5 file, organized by sequence accession IDs and embedding types.

3. **GO Term Lookup**:
   - Embeddings are compared for similarity using a **vector database**.
   - Retrieves GO terms associated with the most similar proteins.
   - Results include GO terms, distances, and metadata.

4. **Results**:
   - Annotations are saved in CSV files for each execution, named with a timestamp for reproducibility.

This pipeline results from joint efforts with equal contribution between **Ana Roja's lab**
(Andalusian Center for Developmental Biology, CSIC) and **Rosa Fernández's lab**
(Metazoa Phylogenomics Lab, Institute of Evolutionary Biology, CSIC-UPF).
This collaboration highlights how synergistic efforts between labs with different expertise can achieve great outcomes.
We thank **LifeHUB-CSIC** for being the catalyst of this project and inspiring us to 'think big'.

Usage
-----
`Setup Instructions <../deployment/setup_instructions.html>`_


Configuration`
---------------------------------------------------------------
`FANTASIA Parameters <../configuration/fantasia/parameters.html>`_


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
import argparse
from datetime import datetime
from protein_metamorphisms_is.sql.model import (
    SequenceEmbedder,
    EmbeddingLookUp,
    AccessionManager,
    UniProtExtractor,
    SequenceEmbeddingManager,
)

from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config


def main(
        config_path="./pipelines/fantasia/config.yaml"):
    conf = read_yaml_config(config_path)
    # AccessionManager(conf).load_accessions_from_csv()
    AccessionManager(conf).fetch_accessions_from_api()
    UniProtExtractor(conf).start()
    SequenceEmbeddingManager(conf).start()
    current_date = datetime.now().strftime("%Y%m%d%H%M%S")

    embedder = SequenceEmbedder(conf, current_date)
    embedder.start()

    lookup = EmbeddingLookUp(conf, current_date)
    lookup.start()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the pipeline with a specified configuration file.")
    parser.add_argument("--config", type=str, required=True, help="Path to the configuration YAML file.")
    args = parser.parse_args()
    main(args.config)

