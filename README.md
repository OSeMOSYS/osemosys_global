# OSeMOSYS global

OSeMOSYS Global is an open-source, open-data model generator for creating 
global energy system models. It can be used to create inter-connected energy 
systems models for both the entire globe and for any geographically diverse 
subset of the globe. Compared to other existing global models, OSeMOSYS Global 
creates a full energy system representation, allows for full user flexibility
in determining the modelling detail and geographic scope, and is built using 
the fully open-source [OSeMOSYS](https://osemosys.readthedocs.io/en/latest/) 
energy system model.


## Table of Contents

- [Documentation](#documentation)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Getting Started](#getting-started)
- [Contributing](#contributing)
- [Support + Feedback](#support--feedback)
- [License](#license)
- [Citing](#citing)

## Documentation
TBD

## Dependencies
OSeMOSYS Global relies on [Snakemake](https://snakemake.readthedocs.io/en/stable/), 
a [Python](https://www.python.org/downloads/) based workflow management system 
to execute the workflow. This requires [Conda](https://docs.conda.io/projects/conda/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html), an open-source package management system. 

OSeMOSYS data is stored using a [Tabular Datapackage](https://specs.frictionlessdata.io/tabular-data-package/) as demonstrated in the [Simplicity example model](https://github.com/OSeMOSYS/simplicity). The python package [otoole](https://github.com/OSeMOSYS/otoole) is 
used to perform data conversions and generate OSeMOSYS datafiles. 
## Installation
1. [Install miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Clone the git repository

```bash
git clone --recurse-submodules https://github.com/OSeMOSYS/osemosys_global.git 
```

If the repository was cloned without the `--recurse-submodules` flag, run the 
commands

```bash
git submodule init
git submodule update 
```

3. Create a conda environment from the supplied environment file 

```bash
/osemosys_global$ conda env create -f workflow/envs/osemosys-global.yaml
```

4. Activate the new `osemosys-global` environment 

```bash
conda activate osemosys-global
```
## Getting Started 

 - Running the workflow 

TBD

 - Modifying the configuration file 

TBD

 - Resetting to defaults 

TBD

**Option One**

4. In the root folder (`osemosys_global/`) run the command `snakemake -c` to execute the workflow

**Option Two**

4. In the root folder (`osemosys_global/`) run the command `snakemake -c --use-conda` to execute the workflow

## Contributing

We appreciate feedback and contribution to this repo! Please see our [contribution guide](CONTRIBUTING.md) for information on how to conribute. 

## Support + Feedback

 - For asking general usage questions, please use the commuity fourm. 
 - For reporting code and data issues, please use the appropiate issue template

## License

OSeMOSYS Global is liscenced under a MIT liscence. Please see our [LICENSE](LICENSE) doc for more information. 

## Citing 

TBD