# Advanced Functionality 

OSeMOSYS Global exposes numerous ways to customize the model. This page describes what 
these customizations are and how to execute them.

## User Defined Data 

A user may want to use the existing structure of OSeMOSYS Global, but sub in 
their own data. This may be open-data that is different from OSeMOSYS Global's 
compiled data, or proprietary data. This can eaisly be done. 

In this example, we show how to change the `SpecifiedAnnualDemand` data of a 
model generated using OSeMOSYS Global

### 1. Generate Template Data 

Instead of running the full workflow, we will only run the workflow up to 
the generation of the input CSV files through the command: 

```bash
snakemake generate_input_data -c
```

### 2. Modify Parameter Data 

At this point, there will be a folder called `results/<scenario>/data/` that 
holds all template data for the scenario. Find the csv file `SpecifiedAnnualDemand.csv`
and change the values in the `VALUE` column. This can be done for any parameter, 
not just SpecifiedAnnualDemand. 

:::{warning}
While any parameter can be updated, do NOT change the set definitions. Sets
are identified by being named in files with all capitals (ie. REGION.csv). 
:::

### 3. Run the Remainder of the Workflow 

Restart the workflow, which will build and solve the model, through the command: 

```bash
snakemake -c
```

## Custom Nodes 

A user may wish to introduce custom nodes to improve the spatial resolution of 
a particular region or country. In this example, we show how to break the 
single node country of Indonesia into seven seperate nodes.  

### 1. Remove the Existing Node

In the configuration file, remove the existing Indonesia node
(`INDXX`), through the `nodes_to_remove` parameter

```yaml
nodes_to_remove:
  - 'IDNXX'
```

### 2. Add New Nodes

In the configuration file, add the names of the new nodes you
wish to add. In this case we will add the nodes `IDNSM`, `IDNJW`, `IDNNU`, 
`IDNKA`, `IDNSL`, `IDNML`, `IDNPP`. 

:::{warning}
The custom nodes being added must start with the 3-letter code of the country 
that they are in (`IDN` in this case). Only the 2-letter regional code can change. 
:::

```yaml
nodes_to_add:
  - 'IDNSM'
  - 'IDNJW'
  - 'IDNNU'
  - 'IDNKA'
  - 'IDNSL'
  - 'IDNML'
  - 'IDNPP'
```

### 3. Add Residual Capacity 

In the folder `resources/data/custom_nodes`, add any residual capacity to the 
`residual_capacity.csv`. The data must be formatted as show. The 
fuel type must be one of the existing fuels in the model (see 
[here](./model-structure.md#acronyms)) and the capacity is given in `GW`.

| FUEL_TYPE | CUSTOM_NODE | START_YEAR | END_YEAR | CAPACITY |
|-----------|-------------|------------|----------|----------|
| COA       | IDNJW       | 2010       | 2040     | 5000     |

### 4. Add Specified Annual Demand

In the folder `resources/data/custom_nodes`, add the annual demand for each 
node to the `specified_annual_demand.csv`. The data must be formatted as 
shown, with the demand given in PJ. 

| CUSTOM_NODE | YEAR | VALUE  |
|-------------|------|--------|
| IDNSM       | 2020 | 319.92 |
| IDNSM       | 2021 | 333.22 |
| IDNSM       | 2022 | 346.53 |

### 5. Add Specified Demand Profile 

In the folder `resources/data/custom_nodes`, add the demand profile for each 
node to the `specified_demand_profile.csv`. The data must be formatted as shown.

| Month | Day | Hour | IDNSM | IDNJW | IDNKA | IDNNU | IDNSL | IDNML | IDNPP |
|-------|-----|------|-------|-------|-------|-------|-------|-------|-------|
| 1     | 1   | 0    | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  |
| 1     | 1   | 1    | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  |
| 1     | 1   | 2    | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  |
| 1     | 1   | 3    | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  |
| 1     | 1   | 4    | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  | 0.04  |

:::{note}
Demand profile data must be given for all hours in a single year. The data 
will then be aggregated according to the timeslice definition and replicated 
for each year. 
:::

### 6. Run the model 

Check that the country with custom nodes is listed under the geographic scope
in the configuration file 

```yaml
geographic_scope:
    - 'IDN'
```

Then run the workflow as normal from the root directory

```bash
snakemake -c
```