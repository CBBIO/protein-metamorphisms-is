Run FANTASIA
============

This guide explains how to run FANTASIA (Functional ANnoTAtion based on embedding space SImilArity) from the command line.

Prerequisites
-------------

Before running FANTASIA, ensure that:

- The software is properly installed.
- All dependencies are met.
- You have a valid configuration file (e.g., `config/constants.yaml`).
- If using the default project structure, the working directory should be set to the root of the project.

Basic Usage
-----------

FANTASIA can be run using the Python module execution syntax. The basic command is:

.. code-block:: bash

   python -m protein_metamorphisms_is.pipelines.fantasia.main --config CONFIG_PATH --input FANTASIA_INPUT_FASTA

Where:

- ``--config CONFIG_PATH`` is the path to your configuration YAML file. If not provided, the default path `./pipelines/fantasia/config/constants.yaml` will be used.
- ``--input FANTASIA_INPUT_FASTA`` is the path to the input FASTA file containing the protein sequences for annotation.

Example
-------

To run FANTASIA with your input FASTA file:

.. code-block:: bash

   python -m protein_metamorphisms_is.pipelines.fantasia.main --input ~/fantasia/input/input.fasta

Option 1: Using Default Project Structure
-----------------------------------------

In this option, you clone the repository, set the working directory to the root of the project, and run FANTASIA with the default configuration and constants files. This is the recommended approach if you want to avoid manually specifying paths.

Steps:

1. Clone the repository and navigate to the root directory:

   .. code-block:: bash

      git clone https://github.com/CBBIO/protein-metamorphisms-is.git
      cd protein-metamorphisms-is

2. Run the pipeline:

   .. code-block:: bash

      python -m protein_metamorphisms_is.pipelines.fantasia.main

3. Place your input FASTA file in the default input directory (`~/fantasia/input/example.fasta`) if you are using this setup.

.. important::

   Ensure that the data directory (by default `~/fantasia`) exists and that you have read and write permissions for this directory.

Option 2: Specifying Custom Paths
---------------------------------

If you prefer more flexibility or do not want to adhere to the default project structure, you can manually specify the paths to the configuration file, constants file, and input FASTA file. This option is useful if you have a custom setup or need to run the pipeline in a different environment.

Steps:

1. Provide the paths to the necessary files explicitly using command-line arguments:

   .. code-block:: bash

      python -m protein_metamorphisms_is.pipelines.fantasia.main --config /path/to/config.yaml --constants /path/to/constants.yaml --input /path/to/input.fasta

Where:

- ``--config /path/to/config.yaml`` is the path to your custom configuration YAML file.
- ``--constants /path/to/constants.yaml`` is the path to your custom constants YAML file.
- ``--input /path/to/input.fasta`` is the path to the input FASTA file containing the protein sequences for annotation.

.. important::

   Make sure that the paths provided are correct and that you have the necessary read and write permissions for the directories.

Pipeline Steps
--------------

When executed, FANTASIA will:

1. Load and validate the configuration.
2. Fetch protein accessions from the API.
3. Extract UniProt data.
4. Generate sequence embeddings.
5. Process the embeddings.
6. Perform GO term lookup.

Output
------

The pipeline generates several output files:

- Embedding files in HDF5 format.
- Results files with timestamps to ensure reproducibility.
- Log files documenting the execution process.

Timestamps are added to output files to prevent overwriting previous results and ensure reproducibility.

Troubleshooting
---------------

If you encounter issues, please:

- Ensure the working directory is set to the root of the project (if using Option 1).
- Check that your configuration and constants files exist and are properly formatted (if using Option 2).
- Verify that all paths in the configuration file are correct.
- Ensure you have sufficient disk space for the output files.
- Check that the data directory (default `~/fantasia`) exists and has proper read/write permissions.
- Check the log files for detailed error messages.

For additional support, please contact the development team using the contact information provided in the main documentation.
