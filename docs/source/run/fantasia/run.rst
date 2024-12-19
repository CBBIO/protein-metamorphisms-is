Run FANTASIA
===========

This guide explains how to run FANTASIA (Functional ANnoTAtion based on embedding space SImilArity) from the command line.

Prerequisites
------------

Before running FANTASIA, ensure that:

1. The software is properly installed
2. All dependencies are met
3. You have a valid configuration file (``config.yaml``)

Basic Usage
----------

FANTASIA can be run using the Python module execution syntax. The basic command is::

    python -m protein_metamorphisms_is.pipelines.fantasia.main [--config CONFIG_PATH]

Parameters
----------

--config CONFIG_PATH
    Path to your configuration YAML file. If not provided, the default path ``./pipelines/fantasia/config.yaml`` will be used.

Example
-------

To run FANTASIA with a specific configuration file::

    python -m protein_metamorphisms_is.pipelines.fantasia.main --config /home/bioxaxi/PycharmProjects/protein-metamorphisms-is/protein_metamorphisms_is/config/config.yaml


Pipeline Steps
-------------

When executed, FANTASIA will:

1. Load and validate the configuration
2. Fetch protein accessions from the API
3. Extract UniProt data
4. Generate sequence embeddings
5. Process the embeddings
6. Perform GO term lookup

Output
------

The pipeline generates several output files:

- Embedding files in HDF5 format
- Results files with timestamps
- Log files documenting the execution process

Timestamps are added to output files to ensure reproducibility and prevent overwriting previous results.

Troubleshooting
--------------

If you encounter issues:

1. Check that your configuration file exists and is properly formatted
2. Verify that all paths in the configuration file are correct
3. Ensure you have sufficient disk space for the output files
4. Check the log files for detailed error messages

For additional support, please contact the development team using the contact information provided in the main documentation.