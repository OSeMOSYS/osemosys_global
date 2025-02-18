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
[Snakemake](https://snakemake.readthedocs.io/en/stable/). A high level overview 
of OSeMOSYS Global's workflow is shown below. The green boxes highlight where 
the user interfaces with the model, while the blue boxes highlight automated 
actions that run behind the scenes.

![Flowchart-high-level](_static/flowchart-high-level.jpg "Flowchart")

The main components of the directory the user will interface with are
highlighted below. This directory structure follows the recommended 
[snakemake structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).
A full overview of the project directory is available in our 
[contributing guidelines](contributing.md#directory-structure).

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
the model and datasets. These configuration options are stored in the 
(conveniently named!)
[configuration file](https://github.com/OSeMOSYS/osemosys_global/tree/master/config). 
An overview of the configuration file options are shown below:

```{include} ../config/README.md
```

:::{seealso}
Our [model structure](./model-structure.md) document for full details on 
technology and commodity names
:::

## Feedback

If you are experiencing issues running any of the examples, please submit a 
[new issue](https://github.com/OSeMOSYS/osemosys_global/issues/new/choose). 
Our GitHub 
[discussion forum](https://github.com/OSeMOSYS/osemosys_global/discussions) is 
also a great place to ask general OSeMOSYS Global questions.
