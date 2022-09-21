# Installation

OSeMOSYS Global is built ontop of the [OSeMOSYS](http://www.osemosys.org/) framework, a linear cost optimization frameowrk. This page will walk through all required installation steps to run an OSeMOSYS Global model from start to finish. 

:::{important}
OSeMOSYS Global is currently only tested on **Linux** systems. If you are using Windows, we suggest you install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install)
:::

## Install GLPK

The GNU GLPK package is used to create the linear programming file. Full installion instructions can be found on their [GLPK](https://www.gnu.org/software/glpk/).

Once installed run the command `glpsol` in the command line. The message below will display indicating that GLPK has installed correctly. 

    ```bash
     ~/osemosys_global$ glpsol

    GLPSOL: GLPK LP/MIP Solver, v4.65
    No input problem file specified; try glpsol --help
    ```

## Install a Solver

OSeMOSYS Global supports three solvers. Included is the open-source solver [CBC](https://github.com/coin-or/Cbc), and the commercial solvers [Gurobi](https://www.gurobi.com/) and [CPLEX](https://www.ibm.com/analytics/cplex-optimizer). **You must install at least one of these solvers**. If you are 
uncertain about which one, we suggest [CBC](https://github.com/coin-or/Cbc) as it is free and the easiest to setup. 

:::{note}
Support for the new open-source solver, [HiGHS](https://highs.dev/), is planned for a later release. 
:::

To verify that **CBC** has installed correctly, type `cbc` into the command line and the following message should be displayed. 

    ```bash
     ~/osemosys_global$ cbc

    Welcome to the CBC MILP Solver 
    Version: 2.10.3 
    Build Date: Mar 24 2020 

    CoinSolver takes input from arguments ( - switches to stdin)
    Enter ? for list of commands or help
    Coin:
    ``` 

To verify that **CPLEX** has installed correctly, type `cplex` into the command line and the following message should be displayed. Type `q` to exit CPLEX.  

    ```bash
     ~/osemosys_global$ cplex

    Welcome to IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer 22.1.0.0
    with Simplex, Mixed Integer & Barrier Optimizers
    5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
    Copyright IBM Corp. 1988, 2022.  All Rights Reserved.

    Type 'help' for a list of available commands.
    Type 'help' followed by a command name for more
    information on commands.

    CPLEX> 
    ``` 

To verify that **Gurobi** has installed correctly, type `gurobi` into the command line and the following message should be displayed.  


## Install Miniconda

OSeMOSYS Global processes data using a series of Python scripts. [Anaconda](https://www.anaconda.com/) is a package manager that manages project dependencies. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is a minimal install of Anaconda. Install miniconda from their [website](https://docs.conda.io/en/latest/miniconda.html)

:::{note}
You may need to restart your terminal for conda to finish installing
:::

To verify that conda is installed, run the command `conda info`. Information about your conda environment will be printed out. 

    ```bash
    (base) ~/osemosys_global$ conda info
    ```

## Clone Repository

All source code for OSeMOSYS Global is hosted on it's own [GitHub repository](https://github.com/OSeMOSYS/osemosys_global), nested under the [OSeMOSYS repository](https://github.com/OSeMOSYS)

    :::{info}
    If you are new to [Git](https://git-scm.com) and [GitHub](https://github.com/), do not worry! All relevent commands are provided below. If you would like to learn more about version control, GitHub has great [documentation](https://docs.github.com/en/get-started/quickstart/hello-world).
    :::

Clone the repository, being sure to include the `--recurse-submodules` flag

    ```bash
     ~/osemosys_global$ git clone --recurse-submodules https://github.com/OSeMOSYS/osemosys_global.git 
    ```

If the repository was cloned without the `--recurse-submodules` flag, run the commands `submodule init` and `submodule update `

    ```bash
    ~/osemosys_global$ git submodule init
    ~/osemosys_global$ git submodule update 
    ```

## Create the Conda Environment

Create a conda environment from the supplied environment file

    ```bash
    (base) ~/osemosys_global$ conda env create -f workflow/envs/osemosys-global.yaml    
    ```

    :::{info}
    The installation of the `osemosys-global` environment may take a few minutes. This is normal.
    :::

Once installed, activate the new `osemosys-global` environment. You should then see `(osemosys-global)` at the start of your command prompt.

    ```bash
    (base) ~/osemosys_global$ conda activate osemosys-global
    ```

    ```bash
    (osemosys-global) ~/osemosys_global$ 
    ```

# Troubleshooting

Sometimes installation doesn't always go as planned... If you are experiencing issues please submit a [new issue](https://github.com/OSeMOSYS/osemosys_global/issues/new/choose) or ask a question on our [discussion fourm](https://github.com/OSeMOSYS/osemosys_global/discussions). 
