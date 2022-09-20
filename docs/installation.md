# Installation

OSeMOSYS Global is built ontop of the [OSeMOSYS](http://www.osemosys.org/) framework, a linear cost optimization frameowrk. This page will walk through all required installation steps to run an OSeMOSYS Global model from start to finish. 

:::{info}
OSeMOSYS Global is currently only tested on Linux systems. If you are using Windows, we suggest you install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install)
:::

## Install GLPK

The GNU GLPK package is used to create the linear programming file. Full installion instructions can be found on their website. 

1. Install [GLPK](https://www.gnu.org/software/glpk/)

## Install a Solver

OSeMOSYS Global supports three solvers. Included is the open-source solver [CBC](https://github.com/coin-or/Cbc), and the commercial solvers [Gurobi](https://www.gurobi.com/) and [CPLEX](https://www.ibm.com/analytics/cplex-optimizer)

:::{info}
Support for the new open-source solver, [HiGHS](https://highs.dev/), is planned for a later release. 
:::

2. Install at least one solver:
a. [CBC](https://github.com/coin-or/Cbc) (open-source)
b. [Gurobi](https://www.gurobi.com/)
c. [CPLEX](https://www.ibm.com/analytics/cplex-optimizer)

## Install Miniconda

OSeMOSYS Global processes data using a series of Python scripts. [Anaconda](https://www.anaconda.com/) is a package manager that manages project dependencies. Miniconda is a minimal install of Anaconda. 

3. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)

## Clone Repository

4. Clone the git repository.

   ```console
    ~/osemosys_global$ git clone --recurse-submodules https://github.com/OSeMOSYS/osemosys_global.git 
    ```

    If the repository was cloned without the `--recurse-submodules` flag, run
    the commands

    ```bash
    ~/osemosys_global$ git submodule init
    ~/osemosys_global$ git submodule update 
    ```

## Create the Conda Environment

5. Create a conda environment from the supplied environment file

    ```bash
    (base) ~/osemosys_global$ conda env create -f workflow/envs/osemosys-global.yaml    
    ```

6. Activate the new `osemosys-global` environment

    ```bash
    (base) ~/osemosys_global$ conda activate osemosys-global
    ```

    ```bash
    (osemosys-global) ~/osemosys_global$ 
    ```

# Troubleshooting

Sometimes installation doesn't always go as planned... If you are experiencing issues with the install, check that GLPK, CBC, and Miniconda were install correctly. 

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