# Installation

OSeMOSYS Global is built ontop of [OSeMOSYS](http://www.osemosys.org/), a linear cost optimization frameowrk. This page will walk through all required installation steps to run an OSeMOSYS Global model from start to finish. 

:::{warning}
OSeMOSYS Global is currently tested on **UNIX** systems. If you are using  Windows, we suggest you install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install)
:::

## Installation Steps

:::{seealso}
If you are new to [Git](https://git-scm.com) and [GitHub](https://github.com/), not to worry! All relevant Git commands are provided below. However, if you would like to learn more about version control, GitHub has great [documentation](https://docs.github.com/en/get-started/quickstart/hello-world).
:::

### 1. Clone Repository

Users can clone the the OSeMOSYS Global [GitHub repository](https://github.com/OSeMOSYS/osemosys_global) using `HTTPS`, `SSH`, or via the GitHub CLI. See [GitHub docs](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) for information on the different cloning methods. 

#### Using `HTTPS`

```bash
$ git clone https://github.com/OSeMOSYS/osemosys_global.git
```

#### Using `SSH`

:::{note}
You will first need to set up a `SSH` key following [instructions by GitHub](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account) to clone using `SSH`. 
:::

```bash 
git clone git@github.com:OSeMOSYS/osemosys_global.git
```

#### Using GitHub CLI

```bash 
gh repo clone OSeMOSYS/osemosys_global
```

### 2. Install Mamba

OSeMOSYS Global uses `conda`/`mamba` to manage project dependencies. We recommend downloading [mamba](https://mamba.readthedocs.io/). Installation instructions to install `mamba` through `Miniforge` can be found on the `mamba` website [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). To verify that mamba is installed, run the command `mamba info`. Information about your mamba environment will be printed out. 

```bash
(base) $ mamba info

     active environment : base
    active env location : /home/xxx/miniforge3
            shell level : 3
       user config file : /home/xxx/.condarc
 populated config files : /home/xxx/miniforge3/.condarc
          conda version : 24.3.0
    conda-build version : not installed
         python version : 3.10.14.final.0
...
```

:::{note}
You may need to restart your terminal for `mamba` to finish installing. Once `mamba` is installed, you will see `(base)` in front at the start of your command line.
:::

### 3. Create the Environment

OSeMOSYS Global stores all project dependencies in a file that mamba can read to create a new environment. Run the command below from the root OSeMOSYS Global directory to create a new environment called `osemosys-global`.

```bash
(base) $ mamba env create -f workflow/envs/osemosys-global.yaml    
```

Once installed, activate the new `osemosys-global` environment. You will now see `(osemosys-global)` at the start of your command prompt.

```bash
(base) $ mamba activate osemosys-global
(osemosys-global) $ 
```

### 4. Install a Solver

OSeMOSYS Global supports three solvers; [`CBC`](https://github.com/coin-or/Cbc), [`Gurobi`](https://www.gurobi.com/) and [`CPLEX`](https://www.ibm.com/analytics/cplex-optimizer). Moreover, OSeMOSYS uses [`GLPK`](https://www.gnu.org/software/glpk/) to generate solver independent linear programming file. To run OSeMOSYS Global, you **must** install `GLPK` and at least one solver. If you installed the dependencies through the environment file, `GLPK` and `CBC` should be automatically installed.  

#### 4.1. Install GLPK

[GNU GLPK](https://www.gnu.org/software/glpk/) is an open-source linear programming package that **will be installed with the environment**. OSeMOSYS Global uses it to create a linear programming files. To confirm that `GLPK` installed correctly, run the command `glpsol` in the command line. The following message should display. 

```bash 
(osemosys-global) $ glpsol

GLPSOL: GLPK LP/MIP Solver, v4.65
No input problem file specified; try glpsol --help
```

:::{note}
If for any reason you need to install `GLPK`, you can do so using the command `mamba install glpk`.
:::

#### 4.2. Install CBC

[`CBC`](https://github.com/coin-or/Cbc) is open-source solver that **will be installed with the environment**. To confirm that `CBC` installed correctly, run the command `cbc` in the command line. The following message should display. Type `quit` to exit `CBC`.

```bash
(osemosys-global) $ cbc

Welcome to the CBC MILP Solver 
Version: 2.10.3 
Build Date: Mar 24 2020 

CoinSolver takes input from arguments ( - switches to stdin)
Enter ? for list of commands or help
Coin:
 ``` 

:::{note}
If for any reason you need to install `CBC`, you can do so using the command `mamba install coin-or-cbc`.
:::

#### 4.3. Install CPLEX

If you are an academic researcher or student, you may qualify for the [academic license](https://www.ibm.com/academic/topic/data-science) of IBM's `CPLEX` optimizer. Else, you will need to purchase a [commercial license](https://www.ibm.com/support/pages/downloading-ibm-ilog-cplex-optimization-studio-v1290). Once installed, run the command `cplex` in the command line. The following message will display indicating that CPLEX has installed correctly. Type `quit` to exit CPLEX.

```bash
(osemosys-global) $ cplex

Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer 22.1.0.0
with Simplex, Mixed Integer & Barrier Optimizers
5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
Copyright IBM Corp. 1988, 2022.  All Rights Reserved.

Type 'help' for a list of available commands.
Type 'help' followed by a command name for more
information on commands.

CPLEX> 
``` 

#### 4.4. Gurobi

If you are an academic researcher or student, you may qualify for the [academic license](https://www.gurobi.com/academia/) of Gurobi's optimizer. Else, you will need to purchase a [commercial license](https://www.gurobi.com/products/gurobi-optimizer/). Once installed, run the command `gurobi_cl` in the command line. The following message will display indicating that Gurobi has installed correctly. 

```bash
(osemosys-global) $ gurobi_cl

Set parameter Username
Set parameter LogFile to value "gurobi.log"
Using license file /home/../gurobi.lic

Usage: gurobi_cl [--command]* [param=value]* filename
Type 'gurobi_cl --help' for more information.
```

### 5. Snakemake 

OSeMOSYS Global uses the workflow management tool [Snakemake](https://snakemake.readthedocs.io/en/stable/) to manage model runs. This should have been installed with the environment. To confirm it has been installed correctly, run the following command from the command line. 

```bash
(osemosys-global) $ snakemake --version 
8.25.3
```

### 6. Run a Model

And thats it! You can now follow [our examples](examples.md) to create an OSeMOSYS model for yourself. 

## Troubleshooting

Sometimes installation doesn't always go as planned... If you are experiencing issues please search our [issue tracker](https://github.com/OSeMOSYS/osemosys_global/issues), and submit an [issue](https://github.com/OSeMOSYS/osemosys_global/issues/new/choose) if the problem is new.

## Dependencies

OSeMOSYS Global relies on numerous open-source community supported tools. Below is a list on the heavily used programs and packages we hope you will investigate further!

[Python](https://www.python.org/downloads/)
: All data processing is writen in the Python programming language

[Mamba](https://mamba.readthedocs.io/en/) 
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
