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

**Option One**

1. [Install miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Create a conda environment with `conda env create -f workflow/envs/osemosys-global.yaml`
3. Activate the new `osemosys-global` environment with `conda activate osemosys-global`
4. In the root folder (`osemosys_global/`) run the command `snakemake -c` to execute the workflow

**Option Two**

1. [Install miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. [Install snakemake](https://snakemake.readthedocs.io/en/stable/#) in a new environment with `conda create -c conda-forge -c bioconda -n osemosys-global snakemake`
3. Activate the new `osemosys-global` environment with `conda activate osemosys-global`
4. In the root folder (`osemosys_global/`) run the command `snakemake -c --use-conda` to execute the workflow
