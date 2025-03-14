from protein_metamorphisms_is.sql.model.entities.embedding.structure_3di import Structure3Di
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.sql.model.entities.structure.chain import Chain
from protein_metamorphisms_is.sql.model.entities.structure.state import State
from protein_metamorphisms_is.sql.model.operational.clustering.cluster import ClusterEntry, Cluster, Subcluster, \
    SubclusterEntry
from protein_metamorphisms_is.tasks.base import BaseTaskInitializer
from protein_metamorphisms_is.helpers.clustering.cdhit import calculate_cdhit_word_length

from pycdhit import cd_hit, read_clstr
import os


class StructuralSubClustering(BaseTaskInitializer):
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

        query = self.session.query(
            Cluster.id.label('cluster_id'),
            ClusterEntry.id.label('cluster_entry_id'),
            Sequence.sequence.label('sequence'),
            Chain.id.label('chain_id'),
            State.id.label('state_id'),
            State.model_id.label('model_identifier'),
            State.file_path.label('model_file_path'),
            Structure3Di.id.label('structure_3di_id'),
            Structure3Di.embedding.label('structure_3di_embedding')
        ).join(
            Cluster, ClusterEntry.cluster_id == Cluster.id
        ).join(
            Sequence, ClusterEntry.sequence_id == Sequence.id
        ).join(
            Chain, Sequence.id == Chain.sequence_id
        ).join(
            State, Chain.id == State.chain_id
        ).join(
            Structure3Di, State.id == Structure3Di.state_id
        ).filter(
            Cluster.id != None
        )

        results = query.all()

        self.logger.info(f"Loaded {len(results)} embedding entries for subclustering")

        # Group results by cluster_id
        clusters = {}
        for result in results:
            cluster_id = result.cluster_id
            if cluster_id not in clusters:
                clusters[cluster_id] = []
            clusters[cluster_id].append(result)

        self.logger.info(f"Grouped into {len(clusters)} clusters")
        return clusters

    def create_fasta(self, cluster_id, embeddings):
        # Use the configuration file to get the directory
        output_dir = os.path.join(self.conf.get('data_directory', '..'), 'cdhit')
        os.makedirs(output_dir, exist_ok=True)

        fasta_path = os.path.join(output_dir, f"cluster_{cluster_id}_embeddings.fasta")
        self.logger.info(f"Writing embeddings to FASTA file at {fasta_path}")
        with open(fasta_path, "w") as fasta_file:
            for emb in embeddings:
                header = f">{emb.structure_3di_id}"  # Use structure_3di_id instead
                sequence = emb.structure_3di_embedding.replace(' ', '')  # Assuming embeddings are stored as strings
                fasta_file.write(f"{header}\n{sequence}\n")
        return fasta_path

    def process(self, cluster_id, embeddings):
        fasta_path = self.create_fasta(cluster_id, embeddings)
        cluster_results = self.run_cd_hit(fasta_path)
        self.store_entry(cluster_id, cluster_results)

    def run_cd_hit(self, fasta_path):
        cdhit_out_path = fasta_path.replace('.fasta', '')
        self.logger.info(f"Running CD-HIT on {fasta_path}")
        cd_hit(
            i=fasta_path,
            o=cdhit_out_path,
            c=self.conf.get('structural_identity_threshold', 0.65),
            d=0,
            sc=1,
            aL=self.conf.get('structural_alignment_coverage', 0.9),
            M=self.conf.get('memory_usage', 1024),
            T=self.conf.get('max_workers', 4),
            g=self.conf.get('most_representative_search', 1),
            n=calculate_cdhit_word_length(self.conf.get('structural_identity_threshold', 0.65), self.logger),
        )
        self.logger.info(f"Reading CD-HIT output from {cdhit_out_path}.clstr")

        # Debugging output of read_clstr
        cluster_results = read_clstr(f"{cdhit_out_path}.clstr")

        if 'cluster' not in cluster_results.columns:
            cluster_results['cluster'] = 0

        self.logger.info(f"CD-HIT cluster results:\n{cluster_results}")
        return cluster_results

    def store_entry(self, cluster_id, cluster_results):
        self.logger.info(f"Storing subcluster data for cluster {cluster_id}")
        subclusters_dict = {}

        # Calcular los tamaños de cada subcluster
        cluster_counts = cluster_results['cluster'].value_counts().to_dict()

        for _, row in cluster_results.iterrows():
            subcluster_identifier = row.get('identifier')
            subcluster_index = row['cluster']
            is_representative = row.get('is_representative', False)

            if subcluster_index not in subclusters_dict:
                # Obtener el tamaño del subcluster actual
                subcluster_size = cluster_counts.get(subcluster_index, 0)
                description = f"(size: {subcluster_size})"
                subcluster = Subcluster(cluster_id=cluster_id, description=description)
                self.session.add(subcluster)
                self.session.flush()  # Obtener el ID generado para el subcluster
                subclusters_dict[subcluster_index] = subcluster.id
                self.logger.info(f"Created new subcluster with ID: {subclusters_dict[subcluster_index]}")

            subcluster_entry = SubclusterEntry(
                subcluster_id=subclusters_dict[subcluster_index],
                structure_3di_id=subcluster_identifier,
                is_representative=is_representative,
                sequence_length=row['size'],
                identity=row['identity']
            )
            self.session.add(subcluster_entry)

        self.session.commit()
        self.logger.info(f"Subcluster data stored for cluster {cluster_id}")
