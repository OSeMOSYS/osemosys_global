# Getting Started

This page will give you an overview of OSeMOSYS Global's workflow, and walk 
through some simple examples. 

:::{seealso}
Ensure you have followed our 
[installation instructions](installation.md#installation) before running a 
model
:::

## Project Overview

OSeMOSYS Global is an open-source, open-data model generator for creating 
global electricity system models. This includes creating interconnected models
for both the entire globe and for any geographically diverse subset of the 
globe. To manage this data processing we use the workflow management tool, 
[Snakemake](https://snakemake.readthedocs.io/en/stable/), to create a 
configurable workflow. A high level overview of OSeMOSYS Global's workflow is 
shown below. The green boxes highlight where the user interfaces with the model, 
while the blue boxes highlight automated actions that run behind the scenes.

![Flowchart-high-level](_static/flowchart_high_level.png "Flowchart")

The main components of the directory the user will interface with are
highlighted below. This directory structure follows the recommended 
[snakemake structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).
A full overview of the project directory is available in our 
[contribution document](contributing.md#directory-structure).

``` bash
osemosys_global
├── docs                      
├── config           
│   ├── config.yaml  # User configurable setup file           
├── resources                    
├── resutls          # Will appear after running 
├── workflow                         
└── ...
```

## Configuration File

OSeMOSYS Global is a configurable model generator, allowing for full user 
flexibility in determining the time slice structure and geographic scope of 
the model and datasets. These configuration options are stored in a 
(conveniently named)
[configuration file](https://github.com/OSeMOSYS/osemosys_global/tree/master/config). 
An overview of the configuration file options are shown below:

```{include} ../config/README.md
```

:::{seealso}
Our [naming conventions](./naming-conventions.md) document for full detials on 
model nomenclature
:::

## Examples

Below are some simple examples you can follow to understand how OSeMOSYS Global
works. 

:::{caution}
Before running any examples, ensure you first follow our [installation instructions](installation.md#installation) and perform these two steps.

1. Navigate to the root OSeMOSYS Global directory in the command line 

    ```bash
    (base) $ cd ~/osemosys_global 

    (base) ~/osemosys_global$ 
    ```

2. Activate the `osemosys-global` conda environment

    ```bash
    (base) ~/osemosys_global$ conda activate osemosys-global

    (osemosys-global) ~/osemosys_global$
    ```
:::

### Example 1

**Goal**: Run the workflow with default settings. This will produce a model 
of India from 2015 to 2050 with 8 timeslices per year, and solve it using CBC.

:::{warning}
If you installed CPLEX or Gurobi instead of CBC, you must first change this in 
the configuration file at `config/config.yaml`
:::

1. Run the command `snakemake -c`. The time to build and solve the model will
vary depending on your computer, but in general, from start to finish the 
workflow will finish within minutes for this example.

    ```bash
    (osemosys-global) ~/osemosys_global$ snakemake -c
    ```

    :::{tip}
    The `-c` command will instruct Snakemake to use all available computer 
    cores. If your want to restrict this, input a number after the `-c` to 
    specify the number of cores. For example, the command `snakemake -c2` will
    run the workflow using 2 cores. See 
    [snakemake's documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#useful-command-line-arguments) 
    for more information on this.   
    :::

2. Navigate to the newly created `results/` folder. All available automatically 
generated results are summarized below. 

    ``` bash
    osemosys_global             
    ├── resutls            # Will appear after running the workflow
    │   ├── data           # Global CSV OSeMOSYS data         
    │   ├── figs           
    │   │   ├── ...        # Global demand projections
    │   ├── india          # Name of scenario
    │   │   ├── data/      # Scenario input CSV data
    │   │   ├── figures            
    │   │   │   ├── GenerationAnnual.html
    │   │   │   ├── GenerationHourly.html
    │   │   │   ├── TotalCapacityAnnual.html
    │   │   │   ├── TransmissionCapacity2050.jpg
    │   │   │   ├── TransmissionFlow2050.jpg
    │   │   ├── result_summaries    
    │   │   │   ├── Capacities.csv
    │   │   │   ├── Generation_By_Node.csv
    │   │   │   ├── Generation.csv
    │   │   │   ├── Metrics.csv
    │   │   │   ├── TradeFlows.csv
    │   │   ├── results/   # Scenario result CSV data
    │   │   ├── india.txt  # Scenario OSeMOSYS data file                       
    └── ...
    ```

    | File    | Description |
    |---------|-------------|
    | `GenerationAnnual.html` | Plot of system level annual generation by technology |
    | `GenerationHourly.html` | Plot of system level technology generation by timeslice |
    | `TotalCapacityAnnual.html` | Plot of system level annual capacity by technology |
    | `Capacities.csv` | Table of nodal level annual capacity by technology |
    | `Generation_By_Node.csv` | Table of nodal level technology generation by timeslice |
    | `Generation.csv` | Table of system level technology generation by timeslice |
    | `Metrics.csv` | Table of system level cost and emission statistics |
    | `TradeFlows.csv` | Table of nodal level electricity trade by timeslice |
    | `TransmissionCapacityXXXX.jpg` | Transmission capacity plot for last year of model |
    | `TransmissionFlowXXXX.jpg` | Transmission flow plot for last year of model |

3. View sytem level capacity and generation results. 

    :::{caution}
    These results are used to showcase the capabilities of OSeMOSYS Global. The
    actual energy supply mix results may need further analysis, such as removing 
    technology bias though implementing resource limits on nuclear.
    :::

    ![Example-1.1](_static/example_1.1.png "Example-1.1")
    ![Example-1.2](_static/example_1.2.png "Example-1.2")

4. View Demand Projections Sample Results for Asia in the file 
`results/figs/Demand projection Asia.jpg`. Grey Dots Represent Historical Country Level Values for Countries in Asia and the Coloured Dots Show Projected Values.

    ![Example-1.3](_static/example_1.3.png "Example-1.3")

### Example 2

**Goal**: Modify the geographic scope, temporal settings, and emission penalty 
of the model. 

The goal of this scenario will be to change the geographic scope to add 
Bangladesh, Bhutan, and Nepal to the model. Moreover, we will change the model
horizon to be from 2020-2040 and increase the number of timeslices per year 
from 8 to 18. Finally, we will ensure crossborder trade is allowed, set the 
emission penalty to $50/T, and create country level result plots.

:::{tip}
All changes described below occur in the `config/config.yaml` file
:::

1. Change the scenario name to BBIN (**B**angladesh, **B**hutan, **I**ndia, 
and **N**epal)

    ```yaml
    scenario: 'BBIN'
    ```

2. Change the geographic scope to include the mentioned countries.

    ```yaml
    geographic_scope:
    - 'IND'
    - 'BGD'
    - 'BTN'
    - 'NPL'
    ```

3. Change the model horizon to be from 2020 to 2040. Both numbers are inclusive.

    ```yaml
    startYear: 2020
    endYear: 2040
    ```

4. Change the number of day parts to represent three even 8 hour segments per 
day. Note that the start number is inclusive, while the end number is exclusive.

    ```yaml
    dayparts:
    D1: [0, 8]
    D2: [8, 16]
    D3: [16, 24]
    ```

5. Change the number of seasons to represent 6 equally spaced days. Both numbers
are inclusive.  

    ```yaml
    seasons:
    S1: [1, 2]
    S2: [3, 4]
    S3: [5, 6]
    S4: [7, 8]
    S5: [9, 10]
    S6: [11, 12]
    ```

    :::{note}
    A timeslice strucutre of 6 seasons and 3 dayparts will result in a model 
    with 18 timeslices per year; 6 representative days each with 3 timeslices. 
    :::

    :::{seealso}
    The [OSeMOSYS documentation](https://osemosys.readthedocs.io/en/latest/index.html)
    has more information on the OSeMOSYS timeslice parameters.  
    :::

6. Ensure the `crossborderTrade` parameter is set to `True`

    ```yaml
    crossborderTrade: True
    ```

7. Change the emission penalty to 50 $/T

    ```yaml
    emission_penalty: 50 
    ```

8. Run the command `snakemake -c`

    ```bash
    (osemosys-global) ~/osemosys_global$ snakemake -c
    ```

    :::{tip}
    If you run into any issues with the workflow, run the command 
    `snakemake clean -c`. This will delete any auto generated files and 
    bring you back to a clean start. 
    :::

9. Navigate to the `results/` folder to view results from this model run.

    :::{tip}
    If you don't change clean the model results, the previous scenario 
    results are saved as long as you change the scenario name. 
    :::

    Notice how under `figures/`, there is now a folder for each country. By 
    setting the `crossborderTrade` parameter to be true, we tell the 
    workflow to create out both system level and country level plots. 

    ``` bash
    osemosys_global             
    ├── resutls            
    │   ├── data           # Global CSV OSeMOSYS data         
    │   ├── figs           
    │   │   ├── ...        # Global demand projections
    │   ├── india          
    │   ├── BBIN          # Name of scenario
    │   │   ├── data/     # Scenario input CSV data
    │   │   ├── figures   
    │   │   │   ├── BTN
    │   │   │   │   ├── GenerationAnnual.html
    │   │   │   │   ├── TotalCapacityAnnual.html         
    │   │   │   ├── BGD
    │   │   │   │   ├── GenerationAnnual.html
    │   │   │   │   ├── TotalCapacityAnnual.html         
    │   │   │   ├── IND
    │   │   │   │   ├── GenerationAnnual.html
    │   │   │   │   ├── TotalCapacityAnnual.html         
    │   │   │   ├── NPL
    │   │   │   │   ├── GenerationAnnual.html
    │   │   │   │   ├── TotalCapacityAnnual.html         
    │   │   │   ├── GenerationAnnual.html
    │   │   │   ├── GenerationHourly.html
    │   │   │   ├── TotalCapacityAnnual.html
    │   │   │   ├── TransmissionCapacity2040.jpg
    │   │   │   ├── TransmissionFlow2040.jpg
    │   │   ├── result_summaries   # Auto generated result tables
    │   │   ├── results/   # Scenario result CSV data
    │   │   ├── BBIN.txt  # Scenario OSeMOSYS data file                       
    └── ...
    ```

10. View system level 2040 hourly generation results by viewing the file
`results/BBIN/figures/GenerationHourly.html`

    :::{caution}
    These results are used to showcase the capabilities of OSeMOSYS Global. The
    actual energy supply mix results may need further analysis, such as removing 
    technology bias though implementing resource limits on nuclear.
    :::

    ![Example-2](_static/example_2.png "Example-2")

3. View system level metrics for this model run by looking at the file 
`results/BBIN/result_summaries/Metrics.csv`

    :::{caution}
    These results are used to showcase the capabilities of OSeMOSYS Global. The
    actual energy supply mix results may need further analysis.
    :::

    | Metric              | Unit                      | Value |
    |---------------------|---------------------------|-------|
    | Emissions	          | Million tonnes of CO2-eq. | 1221  |
    | RE Share            | %                         | 10    |
    | Total System Cost	  | Billion $                 | 1544  |
    | Cost of electricity | $/MWh                     | 20    |
    | Fossil fuel share   | %                         | 4     |


### Example 3

**Goal**: Rerun the BBIN example with new interconnectors.

The goal of this scenario will be to rerun the BBIN scenario 
([example 2](#example-2)), excpet we will tell the model to install three new
electricity interconnectors. In 2025 we will install a 3GW interconnector 
between India and Nepal. Then in 2030 we will install a 1GW and 750MW
interconnector between India and Bhutan and India and Bangladesh respectively. 

1. Change the scenario name

    ```yaml
    scenario: 'BBIN_Interconnector'
    ```

2. Configure the new interconnectors

    ```yaml
    user_defined_capacity:
      TRNINDNONPLXX: [3, 2025]
      TRNINDNEBTNXX: [1, 2030]
      TRNINDEABGDXX: [0.75, 2030]
    ```

    :::{note}
    The last two letters in the region (ie. `NO`, `NE`, and `EA` for India, or
    `XX` for Nepal, Bhutan and Bangladesh) represent the node in each region. 
    See our [naming conventions](naming-conventions.md#naming-conventions) 
    for more information on this. 
    :::

3. View the trade capacity plot in 2040 by looking at the file 
`results/BBIN_interconnector/figures/TransmissionCapacity2040.jpg`

    <img src="_static/example_3.1.jpg" width="400">

4. View the trade flow plot in 2040 by looking at the file 
`results/BBIN_interconnector/figures/TransmissionFlow2040.jpg`

    <img src="_static/example_3.2.jpg" width="400">


### Example 4

**Goal**: Run a World Example

The goal of this scenario is to run a World scenario from 2015 to 2050 with
8 timeslices, graphing results at a system level only

1. Change the scenario name

    ```yaml
    scenario: 'WORLD'
    ```

2. Delete everything under the geographic scope

    ```yaml
    geographic_scope:
    ```

    :::{caution}
    Do **NOT** delete the `geographic_scope:` keyword
    :::

3. Change the model horizon to be from 2015 to 2050. Both numbers are inclusive.

    ```yaml
    startYear: 2015
    endYear: 2050
    ```

4. Reset the temporal parameters back to defaults.

    ```yaml
    dayparts:
      D1: [0, 12]
      D2: [12, 24]

    seasons:
      S1: [12, 1, 2]
      S2: [3, 4, 5]
      S3: [6, 7, 8]
      S4: [9, 10, 11]
    ```

5. Remove the timeshift to set to UTC time. 

    ```yaml
        timeshift: 0
    ```

5. Set the results to only graph at a system level

    ```yaml
    results_by_country: False
    ```

5. Run the command `snakemake -c` 

    :::{warning}
    This scenario will take multiple hours to run using a commercial solver 
    (Gurobi or CPLEX) on a high performance computer.
    :::

    ```bash
    (osemosys-global) ~/osemosys_global$ snakemake -c
    ```

6. View system level results in the `results/WORLD/figures` folder

## Feedback

If you are experienceing issues running any of the examples, please submit a [new issue](https://github.com/OSeMOSYS/osemosys_global/issues/new/choose) on the GitHub. 

:::{seealso}
Our GitHub [discussion fourm](https://github.com/OSeMOSYS/osemosys_global/discussions) is a great place to ask general OSeMOSYS Global questions. OSeMOSYS' [Google Group](https://groups.google.com/g/osemosys) is a good place to ask questions about the OSeMOSYS framework.   
:::