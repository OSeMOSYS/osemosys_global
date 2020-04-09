# OSeMOSYS global
A global power system model generator for OSeMOSYS

## Getting started

After following the installation instructions, type ``snakemake`` on the command
line to run the workflow.

## Data standards

- OSeMOSYS data is stored using a [Tabular Datapackage](https://specs.frictionlessdata.io/tabular-data-package/)
as demonstrated in the [Simplicity example model](https://github.com/OSeMOSYS/simplicity)
- **otoole** is used to perform data conversions and generate OSeMOSYS datafiles from the pre-processed data

# Installation

1. [Install miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Create a conda environment containing a minimal version of snakemake
   `conda create -c conda-forge -c bioconda -n snakemake snakemake-minimal`
3. Install ***otoole** `pip install otoole`
3. Activate the environment `conda activate snakemake`
