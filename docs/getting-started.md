# Getting Started

This page gives an overview of OSeMOSYS Global's workflow and the configuration options available to the user. It assumes the reader has successfully followed the [installation instructions](installation.md). 

## Project Overview

OSeMOSYS Global is an open-source, open-data model generator for creating user-configurable electricity system models. This includes creating interconnected models for both the entire globe and for any geographically diverse subset of the globe. To manage the data pipeline we use the workflow management tool [`snakemake`](https://snakemake.readthedocs.io/en/stable/). A high-level overview of OSeMOSYS Global's `snakemake` workflow is shown below. The green boxes highlight where the user interfaces with the model, while the blue boxes highlight automated actions that run behind the scenes.

![Flowchart-high-level](_static/flowchart-high-level.png "Flowchart")


## Repository Structure 

Upon cloning OSeMOSYS Global, the directory will look similar to tree structure shown below. Highlighted are the directories where the user is expected to interface with OSeMOSYS Global. Specifically, `config/config.yaml` holds general configuration options, `resources/custom_nodes/` holds data that the user can update to calibrate models or perform more indepth scenario analysis, finally, `results/<scenario>/` will hold all result data from a scenario run. Each of these sectors is described in detail below. 

```bash
osemosys_global
├── docs/                      
├── config/           
│   ├── config.yaml     # General configuration options 
├── resources/       
│   ├── custom/         # User defined data 
│   ├── default/        # Default data (DO NOT CHANGE)
├── results/            # Appears after running 
│   ├── scenario_name/  # Holds model data for your scenario
├── workflow/                         
└── ...
```

## Configuration File

OSeMOSYS Global is a configurable model generator, allowing for full user flexibility in determining the time slice structure and geographic scope of the model and datasets. These configuration options are stored in the [`config/config.yaml` file](https://github.com/OSeMOSYS/osemosys_global/tree/master/config). An overview of the configuration file options are given below.

### Top Level Options
Top level options include parameters that do not directly modify the underlying model data. 

```{eval-rst}
.. csv-table::
   :header-rows: 1
   :widths: 16,16,16,16,16
   :file: config_tables/top_level.csv
```

### Temporal Options
OSeMOSYS Global's temporal structure uses a representative day approach. The minimum resolution is representing each year as a single timeslice. The maximum resolution is representing each month as an averaged 24hr day (ie. 12x24=288 timeslices). The temporal scope is between `2015` and `2100`.

```{eval-rst}
.. csv-table::
   :header-rows: 1
   :widths: 16,16,16,16,16
   :file: config_tables/temporal.csv
```

### Spatial Options
OSeMOSYS Global can natively model 164 countries seperated by 265 nodes. Users can select choose to model the entire global, or any number of countries within the world. Furthermore, functionality exists to add/remove to/from the default spatial resolution. 

```{eval-rst}
.. csv-table::
   :header-rows: 1
   :widths: 16,16,16,16,16
   :file: config_tables/spatial.csv
```

### Generator Options
OSeMOSYS Global will automatically aggregate and ingest geo-located powerplants from global datasets, including associated techno-economic parameters. Options to add to this compiled dataset can eaisly be done. Furthermore, modelling assumptions, such as reserve margins, can be modified through the configuration file. 

```{eval-rst}
.. csv-table::
   :header-rows: 1
   :widths: 20,20,20,20,20,20
   :file: config_tables/generators.csv
```

### Transmission Options
OSeMOSYS Global interfaces with the [Gloabl Transmission Database](https://www.sciencedirect.com/science/article/pii/S2352340924003895?dgcid=rss_sd_all) to ingest exising and planned transmission lines. Users can add to this dataset and modify assumptions applied to transmission lines. 

```{eval-rst}
.. csv-table::
   :header-rows: 1
   :widths: 20,20,20,20,20,20
   :file: config_tables/transmission.csv
```

### Storage Options
OSeMOSYS Global interfaces with the [DOE Gloabl Energy Storage Database](https://gesdb.sandia.gov/) to ingest exising and planned storage units. Users can add to this dataset and modify assumptions applied to storage units. 

```{eval-rst}
.. csv-table::
   :header-rows: 1
   :widths: 20,20,20,20,20,20
   :file: config_tables/storage.csv
```

### Policy Options
A core function of energy planning models is expoloring the impacts of different policy. This section describes the configuration options available to model different policies. 

```{eval-rst}
.. csv-table::
   :header-rows: 1
   :widths: 20,20,20,20,20,20
   :file: config_tables/policies.csv
```

## User Defined Data

Often it may be easier to modify/update data through tabular format. In these instances users can update any of the data found in the `resources/data/custom` folder. Moreover, if custom nodes are defined, you must update this data. More explanation can be found on the [Examples page](./examples.md). Here we describe how to format each user defined data csv file. 

### RE Potentials

Renewable Energy potentials allow you to define maximum installable capacity limits for renewable technologies at a per node level. This can be used, for example, to limit renewable installations due to land restrictions. The following information must be provided: 
  - **Fuel Type**: Energy carrier of the technology to limit.
  - **Node**: Node to apply the limit to. (Note: the header says "CUSTOM_NODE", but this can be applied to any default or custom node.) 
  - **Capacity**: Maximum installable capacity in `GW`.

An example limiting onshore wind (`WON`) capacity to 20GW in Canada (`CAN`) British Columbia (`BC`) is given below. This would be added to the file `resources/data/custom/RE_potentials.csv`.

| FUEL_TYPE | CUSTOM_NODE | CAPACITY |
|-----------|-------------|----------|
| WON       | CANBC       | 20       |

### RE Profiles

Renewable Energy potentials allow you to define capacity factors at a per node level. This can be used, for example, if you have a custom solar profile you wish to ingest into the model.  

For Solar PV (`SPV`), Onshore Wind (`WON`), Offshore Wind (`WOF`), and Concentrated Solar Power (`CSP`), the profiles can be defined at a per-hour level. The following information must be provided: 
- **Datetime**: Hour of year to apply profile to.
- **Node**: Node to apply the profile to.
- **Value**: Rated capacity factor as a percentage.

An example defining onshore wind (`WON`) capacity factors in Canada (`CAN`) British Columbia (`BC`) and Alberta (`AB`) is given below. This would be added to the file `resources/data/custom/RE_profiles_WON.csv`.

| Datetime         | CANBC | CANAB |
|------------------|-------|-------|
| 01/01/2015 0:00  | 21.5  | 32.4  |
| 01/01/2015 0:00  | 22.1  | 32    |
| 01/01/2015 0:00  | 22.5  | 32.6  |
| ...              | ...   | ...   |
| 31/12/2015 22:00 | 11.2  | 35    |
| 31/12/2015 23:00 | 9.8   | 35.8  |

For Hydro (`HYD`), the profiles can be defined at a per-month level. 
- **Month**: Month of the year.
- **Node**: Node to apply the profile to.
- **Value**: Rated availability factor as a percentage.

An example defining Hydro capacity factors in Canada (`CAN`) British Columbia (`BC`) and Alberta (`AB`) is given below. This would be added to the file `resources/data/custom/RE_profiles_HYD.csv`. 

| NAME  | M1   | M2   | M3   | M4   | M5   | M6   | M7   | M8 | M9 | M10  | M11  | M12  |
|-------|------|------|------|------|------|------|------|----|----|------|------|------|
| CANBC | 45.6 | 45.8 | 46.8 | 52.4 | 66.3 | 80   | 80   | 80 | 80 | 56.9 | 46.4 | 45.7 |
| CANAB | 32.6 | 40.2 | 44.9 | 49.2 | 55.3 | 73.2 | 76.7 | 79 | 79 | 60.4 | 51.1 | 45.4 |

### Residual Capacity

Residual capacity allows the user to inject exisitng capascity into the system. This is useful as global datasets may not capture smaller projects at individual nodes. This can be done for generation technologies, storage technologies, and transmission technologies. The following information must be provided:
- **Fuel**: Energy carrier of the residual capacity. (Generation Technologies)
- **Technology**: Technology of the residual capacity. (Storage and Transmission Technologies)
- **Node**: Node to apply the capacity to. (Note: the header says "CUSTOM_NODE", but this can be applied to any default or custom node.) 
- **Year**: Year to add residual capacity to.
- **Capacity**: Residual capacity to inject into the system.

An example of adding an additional `1 GW` of residual Hydro (`HYD`) capacity in Canada (`CAN`) British Columbia (`BC`) in `2023` is given below. This would be added to the file `resources/data/custom/residual_capacity.csv`. 

| CUSTOM_NODE | FUEL_TYPE | START_YEAR | END_YEAR | CAPACITY |
|-------------|-----------|------------|----------|----------|
| CANBC       | HYD       | 2023       | 2023     | 1        |

When adding residual transmission and storage capacity, ensure you follow the pre-defined naming conventions described [here](./model-structure.md#technology-codes). For example, adding an additional `1 GW` transmission line between Canada (`CAN`) British Columbia (`BC`) and Alberta (`AB`) in `2023` is given below. This would be added to the file `resources/data/custom/transmission_build_rates.csv`. 

| TRANSMISSION  | START_YEAR | END_YEAR | MAX_BUILD |
|---------------|------------|----------|-----------|
| TRNCANBCCANAB | 2023       | 2023     | 1         |

### Build Rates

Limiting the amount of capacity that can be built is often needed to control over-expansion of a single technology. Build rates can be applied to generation technologies, storage technologies, and transmission technologies. The following information must be provided:
- **Fuel**: Energy carrier of the residual capacity. (Generation Technologies)
- **Technology**: Technology of the residual capacity. (Storage and Transmission Technologies)
- **Node**: Node to apply the capacity to. (Note: the header says "CUSTOM_NODE", but this can be applied to any default or custom node.) 
- **Start Year**: First year to apply expansion limit to 
- **End Year**: Last year to apply expansion limit to 
- **Capacity**: Capacity expansion limit.

An example of limiting the amount of new Hydro (`HYD`) capacity allowed in Canada (`CAN`) British Columbia (`BC`) to `2GW` between `2025` and `2050` is given below. This would be added to the file `resources/data/custom/residual_capacity.csv`. 

| CUSTOM_NODE | FUEL_TYPE | START_YEAR | END_YEAR | CAPACITY |
|-------------|-----------|------------|----------|----------|
| CANBC       | HYD       | 2025       | 2050     | 2        |

When adding transmission and storage capacity build limits, ensure you follow the pre-defined naming conventions described [here](./model-structure.md#technology-codes). For example, limiting transmission line expansion between Canada (`CAN`) British Columbia (`BC`) and Alberta (`AB`) to `2GW` between `2025` and `2050` is given below. This would be added to the file `resources/data/custom/transmission_build_rates.csv`. 

| TRANSMISSION  | START_YEAR | END_YEAR | MAX_BUILD |
|---------------|------------|----------|-----------|
| TRNCANBCCANAB | 2023       | 2023     | 1         |

### Demand

Both the yearly annual demand demand profile at each node can be modified. The user is free to update either or both of these parameters. THis can be useful for calibrating a model, or adjusting loads based on proprietary demand projections.  

To update the annual demand, the following information is needed: 
- **Node**: Node to apply the demand to. (Note: the header says "CUSTOM_NODE", but this can be applied to any default or custom node.) 
- **Year**: Year to apply the demand to. 
- **Value**: Annual demand in `PJ`. 

An example of updating the annual demand in Canada (`CAN`) British Columbia (`BC`) is given below. This would be added to the file `resources/data/custom/specified_annual_demand.csv`. 

| CUSTOM_NODE | YEAR | VALUE |
|-------------|------|-------|
| CANBC       | 2024 | 10.5  |
| CANBC       | 2024 | 10.8  |
| ...         | ...  | ...   |
| CANBC       | 2049 | 17.9  |
| CANBC       | 2050 | 18.1  |

To update the demand profile, the following information is needed:
- **Month**: Month of datetime to apply profile to.
- **Day**: Day of datetime to apply profile to.
- **Hour**: Hour of datetime to apply profile to.
- **Node**: Node to apply the profile to.
- **Value**: Fractional hourly load profile value. 

An example of updating the annual demand in Canada (`CAN`) British Columbia (`BC`) and Alberta (`AB`) is given below. This would be added to the file `resources/data/custom/specified_demand_profile.csv`. 

| Month | Day | Hour | CANBC   |   CANAB |
|-------|-----|------|---------|---------|
| 1     | 1   | 0    | 0.00009 | 0.00010 |
| 1     | 1   | 1    | 0.00009 | 0.00010 |
| 1     | 1   | 2    | 0.00010 | 0.00011 |
| ...                | ...     | ...     |
| 12    | 31  | 22   | 0.00010 | 0.00011 |
| 12    | 31  | 23   | 0.00009 | 0.00011 |

:::{warning}
The demand profile for each node over the year must sum up to 1.0
:::

### Primary Fuel Limits and Prices

Primary fuels, for example Coal (`COA`), Gas (`Gas`), and Oil (`Oil`), can be mined domestically or purchased from the international markets. Users are free to modify the primary fuel reserves in each country (not node!) and the price to import fuel from the interational markets. 

To update domestic fuel resource limits, the following information is needed: 
- **Fuel**: Primary fuel type 
- **Country**: Country to apply the limit to. Note: this is the 3 letter country code without the 2 letter node code! 
- **Year**: Year to apply limit to. The resource limits applied here are not cululative
- **Value**: Yearly resource limit in `PJ`

An example applying coal and gas resources limits to Canada is given in the example below. This would be added to the file `resources/data/custom/fuel_limmits.csv`. 

| FUEL | COUNTRY | VALUE | YEAR |
|------|---------|-------|------|
| GAS  | CAN     | 2000  | 2024 |
| GAS  | COA     | 1500  | 2024 |
| GAS  | CAN     | 2000  | 2025 |
| GAS  | COA     | 1500  | 2025 |

If domestic fuel reserves are exhausted, fuel can be purchased and imported from the international market. The price to purchase fuels can be customized by the user. To do this, the following information is needed:  
- **Fuel**: Primary fuel type 
- **Country**: Country to apply the limit to. Note: this is the 3 letter country code without the 2 letter node code! 
- **Year**: Year to apply limit to. The resource limits applied here are not cululative
- **Value**: Cost of fuel in `$M/PJ`

:::{note}
Additionally, you can input the cost in other units (for example, `$/MMBTU` and include a conversion ratio as an energy content to convert to `$M/PJ`).
:::

An example assigning gas and coal costs for Canada to import is given below. This would be added to the file `resources/data/custom/fuel_prices.csv`. 

| FUEL | COUNTRY | UNIT  | ENERGY_CONTENT | 2020 | 2025 | 2030 | 2040 | 2050 |
|------|---------|-------|----------------|------|------|------|------|------|
| GAS  | CAN     | $M/PJ | 1              | 7.5  | 7.5  | 7.5  | 7.5  | 7.5  |
| COA  | CAN     | $M/PJ | 1              | 4.68 | 4.68 | 4.68 | 4.68 | 4.68 |
| BIO  | CAN     | $M/PJ | 1              | 8    | 8    | 8    | 8    | 8    |
| URN  | CAN     | $M/PJ | 1              | 12   | 12   | 12   | 12   | 12   |

:::{warning}
You must define the costs over the years `2020`, `2025`, `2030`, `2040`, and `2050` as this matches the forcast years from the World Bank Commodity Markets. Intermediate values are interpoloated. 
:::

### Availability Factors

Global level availability factors for each technology can be set in the file `resources/data/custom/availability_factors.csv`. These file comes autopopulated with defaults. An example of how to set `90%` availability factors for Coal (`COA`), Combined Cycle Gas (`CCG`) and Solar (`SPV`) is given below. 

| technology | value   |
|------------|---------|
| COA        | 0.9     |
| CCG        | 0.9     |
| SPV        | 0.9     |

### Operational Life

Global level operational lives for each technology can be set in the file `resources/data/custom/operational_life.csv`. This file comes autopopulated with defaults. An example of how to set operational lives for Coal (`COA`), Combined Cycle Gas (`CCG`) and Solar (`SPV`) is given below. 

| tech | years |
|------|-------|
| COA  | 30    |
| CCG  | 30    |
| SPV  | 20    |

## Custom Nodes

When working with custom nodes you will need to define the centerpoints of the node locations. This is needed to calculate distances between nodes for transmission line expansion. This can be done in the file `resources/data/custom/centerpoints.csv`. For example, adding a new node to represent Washington State (`WA`) in the United States (`USA`) is shown below. In this case, we position the node location close to Seattle, the major demand center in Washington. 

| region  | lat          | log           |
|---------|--------------|---------------|
| USAWA   | 47.594396616 | -122.35969011 |

In addition to the centerpoints, user defined data for demand, capacity, and renewable profiles (among others) must also be added. See the [examples page](./examples.md) for more information.   

## Units

Unless otherwise specified, all results will follow the units given in the table below.

```{eval-rst}
.. csv-table::
   :header-rows: 1
   :widths: 20,20
   :file: data_tables/units.csv
```

## Data Sources 

OSeMOSYS Global compiles data from an array of datasources. Listed below are the main sources used when running OSeMOSYS Global without any modifications. 

```{eval-rst}
.. csv-table::
   :header-rows: 1
   :widths: 20,20,20
   :file: data_tables/sources.csv
```


<!-- Dashboard is not functional right now. -->

<!-- ## Interactive Dashboard

If a user wishes to explore input and output data, OSeMOSYS Global includes 
an interactive dashboard to do this. This dashboard can only be run 
after a successful model run. 

### 1. Ensure the OSeMOSYS Global scenario has been run 

The dashboard is **not** integreated with the snakemake workflow, therefore, 
it will not automatically look for missing input files. Ensure that the 
workflow has been run and that input data CSVs (`results/<scenario>/data/`)
and result file CSVs (`results/<scenario>/results/`) exist. Moreover, 
ensure the scenario name in the configuration file matches the paths to the
result data. 

### 2. Run the dashboard 

Run the file `dashboard.sh` to start the dashboard in a local host with the 
following command. Either open the hyperlink from the command line or 
copy/paste it into a web browser. 

```bash 
bash dashboard.sh
```

### 3. Explore the dashboard 

The dashboard is to visualize the input and result data from an OSeMOSYS Global 
model run. The dashboard consists of 5 tabs; tab 1 is a options tab, and the 
remaining tabs are different visualization options (described below).

:::{warning}
For large models, the dashboard my be slow to respond. Especially for 
parameters/varibales plotted over timeslices, rather than years; for example
the `ProductionByTechnologyAnnual` plot will respond much quicker than the 
``ProductionByTechnology` plot. 
:::

#### 3.1 Options Tab 

Global plotting options for the dashboard. Of note is the Geographic Scope 
radio buttons. If the option is set to System, the legend axis in graphs will 
be by technology. If the option is set to Country or Region, the legend axis 
will be country or region respectively.

#### 3.2 Geographic Overview Tab

Allows the user to visualize what nodes and transmission lines are (and are 
not) included in the model. The user can hover over each node/line to get the 
corresponding node and line name. 

:::{warning}
Any custom nodes added are not included in this map, as latitudes and longitudes 
are not added with custom nodes. All other tabs will correctly incorporate 
custom node data. 
:::

#### 3.3 Input Data Tab

Plots input data. The plot will take on the values described in the options tab. 
The technology fuel dropdown is dynamic, based on the user parameter selection. This
means if a technology or fuel is not listed in the dropdown, there is no 
associated data with it. The plotting options include Area, Line, Stacked Bar 
and Grouped Bar.

#### 3.4 Result Data Tab

Plots result data. The plot will take on the values described in the options tab. 
The technology fuel dropdown is dynamic, based on the user variable selection. This
means if a technology or fuel is not listed in the dropdown, there is no 
associated data with it. The plotting options include Area, Line, Stacked Bar 
and Grouped Bar.

#### 3.5 The Transmission Tab 

Plots transmission line result data. The user can select between plotting 
at a system level, or for individual lines. Moreover, the variables include 
poltting options for either directional flow (inports vs. exports) over the 
line, or total magnitdue of flow. The plotting options include Area, Line, 
Stacked Bar and Grouped Bar.

### 4. Quit the dashboard 

To exit the dashboard, close the tab/browser that has the dashboard open, and
press `ctrl+c` or `cmd+c` from the command line. This will stop the local host.  -->