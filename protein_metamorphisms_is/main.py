from protein_metamorphisms_is.helpers.config.yaml import read_yaml_config
from protein_metamorphisms_is.operation.clustering.sequence_clustering import SequenceClustering
from protein_metamorphisms_is.operation.clustering.sequence_structural_embeddings_subclustering import \
    SequenceStructuralEmbeddingsSubClustering
from protein_metamorphisms_is.operation.embedding.sequence_embedding import SequenceEmbeddingManager
from protein_metamorphisms_is.operation.embedding.structure_embedding import StructureEmbeddingManager
from protein_metamorphisms_is.operation.extraction.accessions import AccessionManager
from protein_metamorphisms_is.operation.extraction.pdb import PDBExtractor
from protein_metamorphisms_is.operation.extraction.uniprot import UniProtExtractor


def main(config_path="config/config.yaml"):
    conf = read_yaml_config(config_path)

    AccessionManager(conf).fetch_accessions_from_api()
    UniProtExtractor(conf).start()
    PDBExtractor(conf).start()
    SequenceEmbeddingManager(conf).start()
    StructureEmbeddingManager(conf).start()
    SequenceClustering(conf).start()
    SequenceStructuralEmbeddingsSubClustering(conf).start()



    # GoMetrics(conf).start()db_inserter_callback
    # OpticsClustering(conf).start()
    # GoPrediction(conf).start()
    # GoPredictionMetricsPerProtein(conf).start()


if __name__ == "__main__":
    main()
