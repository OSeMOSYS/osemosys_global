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
to execute the workflow. This requires [Conda](https://docs.conda.io/projects/conda/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html), a open-source package management system. 

OSeMOSYS data is stored using a [Tabular Datapackage](https://specs.frictionlessdata.io/tabular-data-package/) as demonstrated in the [Simplicity example model](https://github.com/OSeMOSYS/simplicity). The python package [otoole](https://github.com/OSeMOSYS/otoole) is 
used to perform data conversions and generate OSeMOSYS datafiles. 

OSeMOSYS Global uses the open-source GNU Linear Programming Kit, [GLPK](https://www.gnu.org/software/glpk/), and the open-source solver solver, [CBC](https://github.com/coin-or/Cbc).

## Installation
1. Install [GLPK](https://www.gnu.org/software/glpk/)
2. Install [CBC](https://github.com/coin-or/Cbc)
3. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)
4. Clone the git repository

```bash
~/osemosys_global$ git clone --recurse-submodules https://github.com/OSeMOSYS/osemosys_global.git 
```

If the repository was cloned without the `--recurse-submodules` flag, run the 
commands

```bash
~/osemosys_global$ git submodule init
~/osemosys_global$ git submodule update 
```

5. Create a conda environment from the supplied environment file 

```bash
(base) ~/osemosys_global$ conda env create -f workflow/envs/osemosys-global.yaml
```

6. Activate the new `osemosys-global` environment 

```bash
(base) ~/osemosys_global$ conda activate osemosys-global
```
### Troubleshooting
1. Verify that GLPK is installed by running the command `glpsol`.
```bash
(osemosys-global) ~/osemosys_global$ glpsol
```
```bash
GLPSOL: GLPK LP/MIP Solver, v4.65
No input problem file specified; try glpsol --help
```
2. Verify that CBC is installed by running the command `cbc`
```bash
(osemosys-global) ~/osemosys_global$ cbc
```
```bash
Welcome to the CBC MILP Solver 
Version: 2.10.3 
Build Date: Mar 24 2020 

CoinSolver takes input from arguments ( - switches to stdin)
Enter ? for list of commands or help
Coin:
```
Quit the solver with the command `quit`

3. Verify that Conda is installed by running the command `conda info`
```bash
(osemosys-global) ~/osemosys_global$ conda info
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
