from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.operations.base.bioinfo_operator import BioinfoOperatorBase
from protein_data_handler.sql.model import PDBChains, Cluster
from pycdhit import cd_hit, read_clstr
import logging


class CDHit(BioinfoOperatorBase):
    """
    Class for processing protein data using the CD-HIT algorithm.

    Extends BioinfoOperatorBase to utilize its database and configuration handling capabilities.
    This class is specifically tailored for clustering protein sequences using the CD-HIT algorithm.

    :param conf: Configuration dictionary containing necessary parameters.
    :type conf: dict
    """

    def __init__(self, conf):
        """
        Initialize the CDHit class.

        Sets up the configuration, initializes the logger, and prepares the database session.

        :param conf: Configuration dictionary passed to the constructor.
        :type conf: dict
        """
        super().__init__(conf, session_required=True)
        self.logger.info("CDHit instance created")

    def start(self):
        """
        Start the sequence clustering process.

        Orchestrates the process of loading protein chains, creating a FASTA file, and clustering sequences.
        Handles exceptions and logs key events.
        """
        try:
            self.logger.info("Starting CD-HIT clustering process")
            chains = self.load_chains()

            self.create_fasta(chains)
            self.cluster()
            self.logger.info("Clustering process completed successfully")

        except Exception as e:
            self.logger.error(f"Error during clustering process: {e}")
            raise

    def load_chains(self):
        """
        Load protein chains from the database.

        Retrieves all PDBChains records from the database.

        :return: List of PDBChains objects.
        :rtype: list
        """
        self.logger.info("Loading protein chains from the database")
        chains = self.session.query(PDBChains).all()
        return chains

    def create_fasta(self, chains):
        """
        Create a FASTA file from protein chains.

        Writes the protein chains to a FASTA formatted file specified in the configuration.

        :param chains: List of PDBChains objects to write to the FASTA file.
        :type chains: list
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
        Perform sequence clustering using CD-HIT.

        Runs the CD-HIT algorithm on the prepared FASTA file, reads the output cluster file, and stores the results.
        """
        fasta_file_path = self.conf.get('fasta_path', './complete.fasta')
        cdhit_out_path = self.conf.get('cdhit_out_path', './out.clstr')
        self.logger.info(f"Running CD-HIT on {fasta_file_path}")
        cd_hit(
            i=fasta_file_path,
            o=cdhit_out_path,
            c=0.7,
            d=0,
            sc=1,
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


