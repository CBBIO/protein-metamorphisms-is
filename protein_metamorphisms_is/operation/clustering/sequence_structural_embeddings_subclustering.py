from sqlalchemy import func, or_, select
from sqlalchemy.orm import aliased
from protein_metamorphisms_is.tasks.base import BaseTaskInitializer
from protein_metamorphisms_is.sql.model import (
    Subcluster, SubclusterEntry, StructureEmbedding, Model, ClusterEntry,
    Cluster, Sequence, PDBReference, PDBChains, Structure
)
from pycdhit import cd_hit, read_clstr
import os


class SequenceStructuralEmbeddingsSubClustering(BaseTaskInitializer):
    def __init__(self, conf):
        super().__init__(conf)

    def start(self):
        self.logger.info("Starting subclustering process...")
        clusters = self.load_sequences_from_clusters()
        for cluster_id, embeddings in clusters.items():
            if len(embeddings) > 1:  # Ensure the cluster has more than one embedding
                self.logger.info(f"Processing subclustering for cluster {cluster_id} with {len(embeddings)} embeddings")
                self.process(cluster_id, embeddings)
            else:
                self.logger.info(f"Skipping cluster {cluster_id} because it only has {len(embeddings)} embedding(s)")

    def load_sequences_from_clusters(self):
        self.logger.info("Loading sequences from clusters...")
        sequence_alias1 = aliased(Sequence)
        sequence_alias2 = aliased(Sequence)

        query = self.session.query(
            Model.id.label('model_id'),
            Model.model_id.label('model_identifier'),
            Model.file_path.label('model_file_path'),
            Structure.id.label('structure_id'),
            Structure.file_path.label('structure_file_path'),
            PDBChains.id.label('pdb_chain_id'),
            PDBChains.chains,
            sequence_alias1.sequence.label('unified_sequence'),
            PDBChains.sequence_id.label('unified_sequence_id'),
            ClusterEntry.id.label('cluster_entry_id'),
            Cluster.id.label('cluster_id'),
            StructureEmbedding.embedding,
            StructureEmbedding.id.label('structure_embedding_id'),
            func.length(StructureEmbedding.embedding).label('embedding_length'),
            func.length(sequence_alias1.sequence).label('unified_sequence_length')
        ).join(
            Structure, Model.structure_id == Structure.id
        ).outerjoin(
            PDBChains, Structure.id == PDBChains.structure_id
        ).outerjoin(
            sequence_alias1, PDBChains.sequence_id == sequence_alias1.id
        ).outerjoin(
            ClusterEntry, sequence_alias1.id == ClusterEntry.sequence_id
        ).outerjoin(
            Cluster, ClusterEntry.cluster_id == Cluster.id
        ).outerjoin(
            StructureEmbedding, Model.id == StructureEmbedding.model_id
        ).filter(
            StructureEmbedding.embedding != None,
            ClusterEntry.id != None,  # Exclude null cluster_entry_id
            StructureEmbedding.embedding != None,  # Exclude null structure_embeddings
            Cluster.id != None  # Exclude null cluster_id
        )

        results = query.all()
        self.logger.info(f"Loaded {len(results)} embedding entries for subclustering")

        # Group results by cluster_id
        clusters = {}
        for result in results:
            if result.cluster_id not in clusters:
                clusters[result.cluster_id] = []
            clusters[result.cluster_id].append(result)

        self.logger.info(f"Grouped into {len(clusters)} clusters")
        return clusters

    def process(self, cluster_id, embeddings):
        fasta_path = self.create_fasta(cluster_id, embeddings)
        cluster_results = self.run_cd_hit(fasta_path)
        self.store_entry(cluster_id, cluster_results)

    def create_fasta(self, cluster_id, embeddings):
        # Use the configuration file to get the directory
        output_dir = os.path.join(self.conf.get('data_directory', '..'), 'cdhit')
        os.makedirs(output_dir, exist_ok=True)

        fasta_path = os.path.join(output_dir, f"cluster_{cluster_id}_embeddings.fasta")
        self.logger.info(f"Writing embeddings to FASTA file at {fasta_path}")
        with open(fasta_path, "w") as fasta_file:
            for emb in embeddings:
                header = f">{emb.structure_embedding_id}"
                sequence = emb.embedding.replace(' ', '')  # Assuming embeddings are stored as strings
                fasta_file.write(f"{header}\n{sequence}\n")
        return fasta_path

    def run_cd_hit(self, fasta_path):
        cdhit_out_path = fasta_path.replace('.fasta', '')
        self.logger.info(f"Running CD-HIT on {fasta_path}")
        print("threshold",self.conf.get('sequence_identity_threshold', 0.91))
        cd_hit(
            i=fasta_path,
            o=cdhit_out_path,
            c=self.conf.get('sequence_identity_threshold', 0.95),
            d=0,
            sc=1,
            aL=self.conf.get('alignment_coverage', 0.9),
            M=self.conf.get('memory_usage', 1024),
            T=self.conf.get('max_workers', 4),
            g=self.conf.get('most_representative_search', 1)
        )
        self.logger.info(f"Reading CD-HIT output from {cdhit_out_path}.clstr")

        # Debugging output of read_clstr
        cluster_results = read_clstr(f"{cdhit_out_path}.clstr")

        self.logger.info(f"CD-HIT cluster results:\n{cluster_results}")
        return cluster_results

    def store_entry(self, cluster_id, cluster_results):
        self.logger.info(f"Storing subcluster data for cluster {cluster_id}")
        subclusters_dict = {}

        for _, row in cluster_results.iterrows():
            subcluster_identifier = row.get('identifier')
            subcluster_index = row.get('cluster')  # El valor generado por CD-HIT para identificar subclusters
            is_representative = row.get('is_representative', False)

            if subcluster_index not in subclusters_dict:
                # Crear un nuevo subcluster si no existe en el diccionario
                subcluster = Subcluster(cluster_id=cluster_id, description=f"Subcluster {subcluster_index}")
                self.session.add(subcluster)
                self.session.flush()  # Obtén el ID generado para el subcluster
                subclusters_dict[subcluster_index] = subcluster.id
                self.logger.info(f"Created new subcluster with ID: {subclusters_dict[subcluster_index]}")

            subcluster_entry = SubclusterEntry(
                subcluster_id=subclusters_dict[subcluster_index],
                structure_embedding_id=subcluster_identifier,  # Usar identifier como structure_embedding_id
                is_representative=is_representative
            )
            self.session.add(subcluster_entry)

        # Hacer commit de todos los subcluster entries al final
        self.session.commit()
        self.logger.info(f"Subcluster data stored for cluster {cluster_id}")


