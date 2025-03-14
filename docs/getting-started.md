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

Often it may be easier to modify/update data through tabular format. In these instances users can update any of the data found in the `resources/data/custom` folder. Moreover, if custom nodes are defined, you must update this data. 

### RE Profiles

### Residual Capacity

### Demand

### Fuel Limits and Prices

### Build Rates

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