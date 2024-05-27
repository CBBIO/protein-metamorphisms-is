from protein_metamorphisms_is.operations.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import PDBChains, Cluster, ClusterEntry, Sequence
from pycdhit import cd_hit, read_clstr


class CDHit(OperatorBase):
    """
    Class for processing protein data using the CD-HIT algorithm, an efficient algorithm for clustering
    and comparing protein or nucleotide sequences.

    Extends the BioinfoOperatorBase to leverage its database and configuration handling capabilities,
    this class specifically focuses on clustering protein sequences. It facilitates the grouping of
    similar sequences, thereby reducing redundancy and improving the efficiency of subsequent analyses.

    Attributes:
        conf (dict): Configuration dictionary containing necessary parameters for CD-HIT and other operations.

    Methods:
        start: Initiates the process of sequence clustering using CD-HIT.
        load_chains: Loads protein chain data from the database for clustering.
        create_fasta: Creates a FASTA file from the protein chain data.
        cluster: Executes the CD-HIT algorithm and processes the output.

    Usage:
        The CDHit class can be utilized in bioinformatics pipelines for sequence analysis, especially where
        sequence redundancy reduction or efficient sequence comparison is required.
    """

    def __init__(self, conf):
        """
        Initialize the CDHit class.

        Sets up the configuration, initializes the logger, and prepares the database session.

        :param conf: Configuration dictionary passed to the constructor.
        :type conf: dict
        """
        super().__init__(conf)
        self.logger.info("CDHit instance created")

    def start(self):
        """
        Start the protein sequence clustering process using CD-HIT.
        """
        try:
            self.logger.info("Starting CD-HIT clustering process for sequences")
            sequences = self.load_sequences()
            self.create_fasta(sequences)
            self.cluster()
            self.logger.info("Clustering process completed successfully")
        except Exception as e:
            self.logger.error(f"Error during clustering process: {e}")
            raise

    def load_sequences(self):
        """
        Retrieve all unique sequences from the database for clustering.

        Fetches all Sequence records from the database. Each sequence will be used for clustering to reduce redundancy.
        """
        self.logger.info("Loading sequences from the database")
        sequences = self.session.query(Sequence.id, Sequence.sequence).all()
        return sequences

    def create_fasta(self, sequences):
        """
        Generate a FASTA file from a list of sequences.

        Args:
            sequences (list): A list of tuples, each containing a sequence ID and its sequence string.
        """
        fasta_path = self.conf.get('fasta_path', './sequences.fasta')
        self.logger.info(f"Writing sequences to FASTA file at {fasta_path}")
        with open(fasta_path, "w") as fasta_file:
            for seq_id, sequence in sequences:
                header = f">{seq_id}"
                fasta_file.write(f"{header}\n{sequence}\n")

    def cluster(self):
        fasta_file_path = self.conf.get('fasta_path', './sequences.fasta')
        cdhit_out_path = self.conf.get('cdhit_out_path', './out.clstr')

        self.logger.info(f"Running CD-HIT on {fasta_file_path}")
        cd_hit(
            i=fasta_file_path,
            o=cdhit_out_path,
            c=self.conf.get('sequence_identity_threshold', 0.9),
            d=0,
            sc=1,
            aL=self.conf.get('alignment_coverage', 0.9),
            M=self.conf.get('memory_usage', 1024),
            T=self.conf.get('max_workers', 4),
            g=self.conf.get('most_representative_search', 1)
        )

        self.logger.info(f"Reading CD-HIT output from {cdhit_out_path}.clstr")
        df_clstr = read_clstr(f"{cdhit_out_path}.clstr")
        clusters_dict = {}
        for _, row in df_clstr.iterrows():
            cluster_id = row['cluster']
            if cluster_id not in clusters_dict:
                cluster = Cluster()
                self.session.add(cluster)
                self.session.flush()  # To obtain the generated cluster ID
                clusters_dict[cluster_id] = cluster.id

            cluster_entry = ClusterEntry(
                cluster_id=clusters_dict[cluster_id],
                sequence_id=row["identifier"],
                is_representative=row['is_representative'],
                sequence_length=row['size'],
                identity=row['identity']
            )
            self.session.add(cluster_entry)
        self.session.commit()
        self.logger.info("CD-HIT clustering data stored in the database")

