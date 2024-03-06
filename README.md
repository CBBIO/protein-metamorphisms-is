[![codecov](https://codecov.io/gh/CBBIO/protein-metamorphisms-is/graph/badge.svg?token=mtOqdG0xbU)](https://codecov.io/gh/CBBIO/protein-metamorphisms-is)
[![PyPI version](https://badge.fury.io/py/protein-metamorphisms-is.svg)](https://badge.fury.io/py/protein-metamorphisms-is)
[![Documentation Status](https://readthedocs.org/projects/protein-metamorphisms-is/badge/?version=latest)](https://protein-metamorphisms-is.readthedocs.io/es/latest/?badge=latest)
![Linting Status](https://github.com/CBBIO/protein-metamorphisms-is/actions/workflows/test-lint.yml/badge.svg?branch=main)


# Information System for Protein Metamorphisms discovery

## Overview
Information System for Protein Metamorphisms discovery is a Python library aimed at simplifying the processing of protein data. It integrates with UniProt and Protein Data Bank (PDB) for searching, downloading, and storing protein information, focusing on metamorphisms discovery in protein data.

## Key Features
- **UniProt and PDB Integration**: Easy access to comprehensive protein data.
- **Database Management**: Supports PostgreSQL for efficient data storage and querying.
- **Flexible Configuration**: User-friendly setup via YAML files.
- **Robust Data Processing**: Utilizes CD-HIT for grouping sequences and supports structural alignment for detailed analysis.

## Getting Started
To begin using the Protein Data Handler, clone the repository and follow the installation instructions in the documentation.

```bash
git clone https://github.com/CBBIO/protein-metamorphisms-is.git
cd protein-data-handler
poetry install
python3 protein_data_handler/main.py
