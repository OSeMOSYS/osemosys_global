# Installation

OSeMOSYS Global is built ontop of [OSeMOSYS](http://www.osemosys.org/), a 
linear cost optimization frameowrk. This page will walk through all required 
installation steps to run an OSeMOSYS Global model from start to finish. 

:::{important}
OSeMOSYS Global is currently tested on **Linux** systems. If you are using 
Windows, we suggest you install 
[Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install)
:::

## Installation Steps

###  1. Install GLPK

The GNU GLPK package is a open-source linear programming package. OSeMOSYS 
Global uses it to create the linear programming file. Full installion 
instructions can be found on their [website](https://www.gnu.org/software/glpk/). 
Once installed run the command `glpsol` in the command line. The following
message will display indicating that GLPK has installed correctly. 

``` bash
$ glpsol

GLPSOL: GLPK LP/MIP Solver, v4.65
No input problem file specified; try glpsol --help
```

### 2. Install a Solver

OSeMOSYS Global supports three solvers. Included is the open-source solver 
[CBC](https://github.com/coin-or/Cbc), and the commercial solvers 
[Gurobi](https://www.gurobi.com/) and 
[CPLEX](https://www.ibm.com/analytics/cplex-optimizer). **You must install at 
least one of these solvers**. If you are uncertain about which one, we suggest 
[CBC](https://github.com/coin-or/Cbc) as it is free and the easiest to setup. 

:::{note}
Support for the new open-source solver, [HiGHS](https://highs.dev/), is planned for a later release. 
:::

#### 2.1. Install CBC

Follow the [download instruction](https://github.com/coin-or/Cbc#download) on 
CBC's GitHub. Once installed, run the command `cbc` in the command line. The 
following message will display indicating that CBC has installed correctly. 
Type `quit` to exit the program.

```bash
$ cbc

Welcome to the CBC MILP Solver 
Version: 2.10.3 
Build Date: Mar 24 2020 

CoinSolver takes input from arguments ( - switches to stdin)
Enter ? for list of commands or help
Coin:
 ``` 

#### 2.2. Install CPLEX

If you are an academnic researcher or student, you may qualify for the 
[academic license](https://www.ibm.com/academic/topic/data-science) of IBM's 
CPLEX optimizer. Else, you will need to purchase a 
[commercial license](https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v1290). 
Once installed, run the command `cplex` in the command line. The following 
message will display indicating that CPLEX has installed correctly. 
Type `quit` to exit the program. Type `q` to exit CPLEX.  

```bash
$ cplex

Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer 22.1.0.0
with Simplex, Mixed Integer & Barrier Optimizers
5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
Copyright IBM Corp. 1988, 2022.  All Rights Reserved.

Type 'help' for a list of available commands.
Type 'help' followed by a command name for more
information on commands.

CPLEX> 
``` 

#### 2.3. Gurobi

If you are an academnic researcher or student, you may qualify for the 
[academic license](https://www.gurobi.com/academia/) of Gurobi's optimizer. 
Else, you will need to purchase a 
[commercial license](https://www.gurobi.com/products/gurobi-optimizer/). 
Once installed, run the command `gurobi_cl` in the command line. The following 
message will display indicating that Gurobi has installed correctly. 

```bash

$ gurobi_cl

Set parameter Username
Set parameter LogFile to value "gurobi.log"
Using license file /home/../gurobi.lic

Usage: gurobi_cl [--command]* [param=value]* filename
Type 'gurobi_cl --help' for more information.
```

### 3. Clone Repository

All source code for OSeMOSYS Global is hosted on it's own 
[GitHub repository](https://github.com/OSeMOSYS/osemosys_global), nested under 
the [OSeMOSYS organization](https://github.com/OSeMOSYS).

:::{seealso}
If you are new to [Git](https://git-scm.com) and [GitHub](https://github.com/), 
not to worry! All relevent Git commands are provided below. However, if you 
would like to learn more about version control, GitHub has great 
[documentation](https://docs.github.com/en/get-started/quickstart/hello-world).
:::

In the command line, navigate to the directory you want to clone the repository 
into. Clone the repository using the command below, being sure to include the
`--recurse-submodules` flag. Navigate into the new `osemosys_global` directory.

```bash
$ git clone --recurse-submodules https://github.com/OSeMOSYS/osemosys_global.git 

$ cd osemosys_global

~/osemosys_global$ 

```

:::{caution}
If the repository was cloned without the `--recurse-submodules` flag, run the 
commands `submodule init` and `submodule update `

```bash
~/osemosys_global$ git submodule init
~/osemosys_global$ git submodule update 
```
:::

### 4. Install Miniconda

OSeMOSYS Global processes and plots data using a series of Python scripts. 
OSeMOSYS Global uses [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
to helps manage all required Python packages. Install Miniconda following their 
[instructions](https://docs.conda.io/en/latest/miniconda.html).

:::{note}
You may need to restart your terminal for conda to finish installing. Once 
Miniconda is istalled, you will see `(base)` infront at the start of your
command line.
:::

To verify that conda is installed, run the command `conda info`. Information 
about your conda environment will be printed out. 

```bash
(base) ~/osemosys_global$ conda info

active environment : base
active env location : ~/miniconda3/envs/base
shell level : 2
user config file : ~/.condarc
populated config files : 
conda version : 4.12.0
...
```

### 5. Create the Conda Environment

OSeMOSYS Global records all project dependencies in a file that conda can read 
from to create a new environment. Run the command below to create a new 
envirnoment called `osemosys-global`. This new envirnoment will install all 
project dependencies. 

```bash
(base) ~/osemosys_global$ conda env create -f workflow/envs/osemosys-global.yaml    
```

:::{caution}
The installation of the `osemosys-global` environment may take a few minutes. This is normal.
:::

Once installed, activate the new `osemosys-global` environment. You will now see 
`(osemosys-global)` at the start of your command prompt.

```bash
(base) ~/osemosys_global$ conda activate osemosys-global

(osemosys-global) ~/osemosys_global$ 
```

### 6. Run a Model

And thats it! You can now follow [our examples](getting-started.md#examples) 
to create a model for yourself. 

## Troubleshooting

Sometimes installation doesn't always go as planned... If you are experiencing 
issues please submit a 
[new issue](https://github.com/OSeMOSYS/osemosys_global/issues/new/choose). 

:::{seealso}
Our GitHub [discussion fourm](https://github.com/OSeMOSYS/osemosys_global/discussions) is a great place to ask general OSeMOSYS Global questions. OSeMOSYS' [Google Group](https://groups.google.com/g/osemosys) is a good place to ask questions about the OSeMOSYS framework.   
:::

## Dependencies

OSeMOSYS Global relies on numerous open-source community supported tools.
Below is a list on the heavily used programs and packages we hope you will 
investigate further!

[Python](https://www.python.org/downloads/)
: All data processing is writen in the Python programming language

[Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
: Python package management tool

[Snakemake](https://snakemake.readthedocs.io/en/stable/)
: Python based workflow management tool

[otoole](https://github.com/OSeMOSYS/otoole)
: Python based command line interface used to work with OSeMOSYS data

[Pandas](https://pandas.pydata.org/) 
: Python package used to transform and analyze data

[Plotly](https://plotly.com/)
: Python package for data visualization

[GLPK](https://www.gnu.org/software/glpk/) 
: Open-source linear programming toolkit

[CBC](https://github.com/coin-or/Cbc)
: Open-source linear program solver
