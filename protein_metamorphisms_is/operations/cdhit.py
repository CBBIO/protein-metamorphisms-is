from protein_metamorphisms_is.operations.base.operator import OperatorBase
from protein_metamorphisms_is.sql.model import PDBChains, Cluster
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

        Coordinates the steps for loading protein chains, creating a FASTA file, executing CD-HIT for clustering, and
        handling exceptions and key events. Logs the progress and any errors encountered during the process.
        """
        try:
            self.logger.info("Starting CD-HIT clustering process")
            chains = self.load_chains()
            # self.explore_representatives()

            self.create_fasta(chains)
            self.cluster()
            self.logger.info("Clustering process completed successfully")

        except Exception as e:
            self.logger.error(f"Error during clustering process: {e}")
            raise

    def load_chains(self):
        """
        Retrieve protein chain data from the database.

        Fetches all PDBChains records from the database. The method can be configured to include or exclude multiple chain
        models based on the 'allow_multiple_chain_models' (NMR samples) configuration.

        Returns:
            list: A list of PDBChains objects representing protein chains.
        """
        self.logger.info("Loading protein chains from the database")
        if not self.conf.get("allow_multiple_chain_models"):
            chains = self.session.query(PDBChains).filter(PDBChains.model == 0).all()
        else:
            chains = self.session.query(PDBChains).all()
        return chains

    def create_fasta(self, chains):
        """
        Generate a FASTA file from a list of protein chains.

        Writes the provided protein chains to a FASTA formatted file. The path for the FASTA file is specified in the
        configuration.

        Args:
            chains (list): A list of PDBChains objects to be written to the FASTA file.
        """
        fasta_path = self.conf.get('fasta_path', './complete.fasta')
        self.logger.info(f"Writing protein chains to FASTA file at {fasta_path}")
        with open(fasta_path, "w") as fasta_file:
            for chain in chains:
                header = f">{chain.id}"
                sequence = chain.sequence
                fasta_file.write(f"{header}\n{sequence}\n")

    def cluster(self):
        """
        Execute the CD-HIT algorithm for sequence clustering.

        Runs the CD-HIT algorithm on the prepared FASTA file, then reads the output cluster file to store the clustering
        results in the database. Configuration parameters such as sequence identity threshold, alignment coverage, accurate mode and
        memory usage are used to control the CD-HIT execution.
        """
        fasta_file_path = self.conf.get('fasta_path', './complete.fasta')
        cdhit_out_path = self.conf.get('cdhit_out_path', './out.clstr')

        sequence_identity_threshold = self.conf.get('sequence_identity_threshold')
        alignment_coverage = self.conf.get('alignment_coverage')
        memory_usage = self.conf.get('memory_usage')
        num_threads = self.conf.get('max_workers')
        most_representative_search = self.conf.get('most_representative_search')

        self.logger.info(f"Running CD-HIT on {fasta_file_path}")
        cd_hit(
            i=fasta_file_path,
            o=cdhit_out_path,
            c=sequence_identity_threshold,
            d=0,
            sc=1,
            aL=alignment_coverage,
            M=memory_usage,
            T=num_threads,
            g=most_representative_search
        )

        self.logger.info(f"Reading CD-HIT output from {cdhit_out_path}.clstr")
        df_clstr = read_clstr(f"{cdhit_out_path}.clstr")
        for _, row in df_clstr.iterrows():
            chain_id = row["identifier"]
            cluster = Cluster(
                cluster_id=row['cluster'],
                pdb_chain_id=chain_id,
                is_representative=row['is_representative'],
                sequence_length=row['size'],
                identity=row['identity']
            )
            self.session.add(cluster)
        self.session.commit()
        self.logger.info("CD-HIT clustering data stored in the database")
