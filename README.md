## Open-Source, Open-Data, Global Electricity System Models

OSeMOSYS Global is an open-source, open-data model generator for creating
global energy system models. It can be used to create inter-connected energy
systems models for both the entire globe and for any geographically diverse
subset of the globe. Compared to other existing global models, OSeMOSYS Global
creates a full energy system representation, allows for full user flexibility
in determining the modelling detail and geographic scope, and is built using
the fully open-source [OSeMOSYS](https://osemosys.readthedocs.io/en/latest/)
energy system model.

## Dependencies

OSeMOSYS Global relies on numerous open-source community supported tools.
Below is a list on the heavily used packages that we hope you investigate
further for yourself!

- [Python](https://www.python.org/downloads/) is used for all data processing
- [Anaconda](https://docs.conda.io/projects/conda/en/latest/) and
[Miniconda](https://docs.conda.io/en/latest/miniconda.html) are used to manage
Python packages
- [Snakemake](https://snakemake.readthedocs.io/en/stable/) is a Python based
workflow management tool
- [otoole](https://github.com/OSeMOSYS/otoole) is a Python based command line
interface used generate OSeMOSYS input and output datafiles
- [Pandas](https://pandas.pydata.org/) is a Python package used to transform
and analyze data
- [Plotly](https://plotly.com/) is a Python package for data visualization
- [GLPK](https://www.gnu.org/software/glpk/) is used to create linear programming
data files
- [CBC](https://github.com/coin-or/Cbc) is a linear programming and mixed integer
program solver
