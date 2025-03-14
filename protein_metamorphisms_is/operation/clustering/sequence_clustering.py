"""
Sequence Clustering Tasks
==========================

The `SequenceClustering` class is responsible for clustering protein sequences using the CD-HIT algorithm.
It extends the `BaseTaskInitializer` class, providing the necessary framework for managing sequence clustering
tasks within the system.

**Purpose**

The `SequenceClustering` class manages the process of clustering similar protein sequences to reduce redundancy
and facilitate further analysis. This is achieved by loading sequences from the database,
running the CD-HIT clustering algorithm, and storing the resulting clusters back into the database.

**Customization**

To create a custom sequence clustering task, subclass `SequenceClustering` and implement any additional methods
or overrides necessary for specific clustering requirements.

**Key Features**

- **Sequence Loading**: Retrieves protein sequences from the database for clustering.
- **FASTA File Generation**: Converts sequences into a FASTA file format required by CD-HIT.
- **CD-HIT Execution**: Runs the CD-HIT algorithm to cluster sequences based on sequence identity thresholds.
- **Cluster Storage**: Stores the resulting clusters and their members into the database.

**Installation of CD-HIT**

Before using the `SequenceClustering` class, the CD-HIT algorithm must be installed on your system.
CD-HIT can be easily installed on Debian-based systems using `apt-get`:

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install cd-hit

This ensures that the CD-HIT executable is available in your system's PATH, which is necessary
for the `SequenceClustering` class to function properly.

**SQL/ORM Entities**

The module interacts with several SQL/ORM entities for organizing and storing clustering results:

- :class:`~protein_metamorphisms_is.sql.model.PDBChains`: Manages Protein Data Bank chains information.
- :class:`~protein_metamorphisms_is.sql.model.Cluster`: Stores clustering results, including cluster identifiers and sequence information.

These entities facilitate the storage and retrieval of clustered protein sequence data in a relational database.

**Configuration**

`SequenceClustering` requires specific configuration parameters for optimal operation. Here's a configuration template:

.. code-block:: yaml

    # Sequence Clustering Configuration
    max_workers: [Number of threads for CD-HIT]
    fasta_path: [Path for FASTA file]
    cdhit_out_path: [Path for CD-HIT output file]
    sequence_identity_threshold: [Identity threshold for clustering]
    memory_usage: [Maximum memory usage for CD-HIT]
    alignment_coverage: [Minimum alignment coverage for clustering]
    most_representative_search: [Boolean value to enable/disable most representative search]

Adjust these settings based on the specific requirements of your clustering tasks.

**Example Usage**

Here is an example of how to subclass `SequenceClustering`:

.. code-block:: python

   from protein_metamorphisms_is.tasks.sequence_clustering import SequenceClustering

   class MySequenceClustering(SequenceClustering):
       def start(self):
           # Custom start logic
           pass
"""
from protein_metamorphisms_is.sql.model.entities.sequence.sequence import Sequence
from protein_metamorphisms_is.sql.model.operational.clustering.cluster import Cluster, ClusterEntry
from protein_metamorphisms_is.tasks.base import BaseTaskInitializer
from pycdhit import cd_hit, read_clstr
from protein_metamorphisms_is.helpers.clustering.cdhit import calculate_cdhit_word_length
import os


class SequenceClustering(BaseTaskInitializer):
    """
    The SequenceClustering class handles the clustering of protein sequences using
    the CD-HIT algorithm. It facilitates the loading of sequences from the database,
    execution of the clustering algorithm, and storage of the clustering results.

    Attributes:
        logger (Logger): Inherited from BaseTaskInitializer for logging purposes.
    """

    def __init__(self, conf):
        """
        Initialize the SequenceClustering.

        This constructor sets up the logger and initializes the configuration.

        Args:
            conf (dict): Configuration dictionary containing settings for clustering.
        """
        super().__init__(conf)
        self.logger.info("CDHit instance created")

    def start(self):
        """
        Start the protein sequence clustering process using CD-HIT.

        This method orchestrates the entire process of loading sequences,
        running the CD-HIT clustering algorithm, and storing the results.
        """
        try:
            self.logger.info("Starting CD-HIT clustering process for sequences")
            sequences = self.load_sequences()

            self.create_fasta(sequences)
            cluster_data = self.process()
            self.store_entry(cluster_data)
            self.logger.info("Clustering process completed successfully")
        except Exception as e:
            self.logger.error(f"Error during clustering process: {e}")
            raise

    def load_sequences(self):
        """
        Load protein sequences from the database.

        This method queries the database to retrieve all protein sequences
        that need to be clustered.

        Returns:
            list: A list of tuples containing sequence IDs and sequences.
        """
        self.logger.info("Loading sequences from the database")
        sequences = self.session.query(Sequence.id, Sequence.sequence).all()
        return sequences

    def create_fasta(self, sequences):
        """
        Create a FASTA file from the loaded sequences.

        This method writes the sequences to a FASTA file, which is the input format
        required by the CD-HIT algorithm.

        Args:
            sequences (list): A list of tuples containing sequence IDs and sequences.
        """
        fasta_path = os.path.expanduser(self.conf.get('fasta_path', './sequences.fasta'))
        self.logger.info(f"Writing sequences to FASTA file at {fasta_path}")
        with open(fasta_path, "w") as fasta_file:
            for seq_id, sequence in sequences:
                header = f">{seq_id}"
                fasta_file.write(f"{header}\n{sequence}\n")

    def process(self):
        """
        Run the CD-HIT algorithm on the sequences.

        This method executes the CD-HIT clustering algorithm on the sequences,
        using the configuration parameters provided.

        Returns:
            DataFrame: A DataFrame containing the clustering results.
        """
        fasta_file_path = os.path.expanduser(self.conf.get('fasta_path', './sequences.fasta'))
        cdhit_out_path = os.path.expanduser(self.conf.get('cdhit_out_path', './out.clstr'))
        self.logger.info(f"Running CD-HIT on {fasta_file_path}")
        cd_hit(
            i=fasta_file_path,
            o=cdhit_out_path,
            c=self.conf.get('sequence_identity_threshold', 0.5),
            d=0,
            sc=1,
            aL=self.conf.get('alignment_coverage', 0.9),
            M=self.conf.get('memory_usage', 1024),
            T=self.conf.get('max_workers', 4),
            g=self.conf.get('most_representative_search', 1),
            n=calculate_cdhit_word_length(self.conf.get('sequence_identity_threshold', 0.5), self.logger),
        )
        self.logger.info(f"Reading CD-HIT output from {cdhit_out_path}.clstr")
        return read_clstr(f"{cdhit_out_path}.clstr")

    def store_entry(self, cluster_data):
        """
        Store the clustering results into the database.

        This method takes the clustering results from CD-HIT and stores each cluster
        and its entries into the database.

        Args:
            cluster_data (DataFrame): A DataFrame containing the clustering results.
        """
        clusters_dict = {}
        for _, row in cluster_data.iterrows():
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
        self.logger.info("Cluster data stored in the database")
