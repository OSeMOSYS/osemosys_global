# OSeMOSYS GNU MathProg

Thanks for using OSeMOSYS and welcome to the OSeMOSYS community.

To use OSeMOSYS you will first need to install a few dependencies. Instructions for
these are given below.

## Running OSeMOSYS

If you've followed the installation instructions below using conda, before running OSeMOSYS you'll
need to activate your conda environment:

    conda activate osemosys

To run OSeMOSYS, first navigate to the folder containing your OSeMOSYS files.
By default, OSeMOSYS results are written into the `results` directory. If this doesn't already exist,
you'll want to create it:

    mkdir results

Now enter the following line into your command prompt:

    glpsol -m osemosys.txt -d  simplicity.txt

The results from the model run should be available in the `results` folder.

### Short and Fast versions of the code

Alternatively, you can run the "short" or "fast" versions of the OSeMOSYS code.
These versions of the OSeMOSYS code produce the same results as the original version, but are written
using a more concise form which condenses many of the equations. For larger models that are more
computationally expensive, the short and fast versions will provide a significant performance boost.
However, for development purposes, the original version of the OSeMOSYS code is much easier to understand. 

    glpsol -m osemosys_short.txt -d  simplicity.txt
    glpsol -m osemosys_fast.txt -d  simplicity.txt

### Running using other solvers

Running using CBC or CLP is a little more involved, as it requires two separate steps.
First, an LP file is generated using glpsol:

    glpsol -m osemosys_short.txt -d simplicity.txt --wlp simplicity.lp --check

Then the LP file is passed to the CBC (or CLP) solver to produce a solution file:

    cbc simplicity.lp solve -solu simplicity.sol

The solution file is written out in a format that requires post-processing. We recommend using
the open-source Python package **otoole** to do this for you:

    otoole results cbc csv simplicity.sol results --input_datafile simplicity.txt

## Installation

To use OSeMOSYS, you'll need to install the 
[GNU Linear Programming Kit](https://www.gnu.org/software/glpk/).

The easiest way to do this is to use Anaconda or Miniconda.

### 1. Install conda package manager

First, install [miniconda](https://docs.conda.io/en/latest/miniconda.html#) by following
the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/#).

### 2. Create a new environment

Run the following commands to create and activate a new environment containing the `glpk` library:

    conda create -n osemosys python=3.8 glpk
    conda activate osemosys

### 3. Optional - install **otoole** for pre- and post- processing

First, activate the conda environment and then install otoole::

    conda activate osemosys
    pip install otoole

### 4. Optional - install CBC or CLP open-source solvers

On Linux or OSX you can also install the more powerful CBC or CLP
open-source solvers using conda:

    conda install -c conda-forge coincbc coin-or-clp

## Where to get help

### The original paper

The [original OSeMOSYS paper](https://doi.org/10.1016/j.enpol.2011.06.033) provides a plain English
description for each of the equation blocks in the model file (`osemosys.txt`).

### OSeMOSYS Documentation

You can read the [OSeMOSYS documentation online](https://osemosys.readthedocs.io/en/latest/?badge=latest).

### Forum

The [Google Groups OSeMOSYS forum](https://groups.google.com/u/1/g/osemosys) 
contains a large number of questions and answers and is useful and active
source of information for how to use OSeMOSYS and tackle all sorts of
OSeMOSYS modelling issues.

### Github Issues

The OSeMOSYS source code and issue tracker are located on [Github](https://github.com/OSeMOSYS/OSeMOSYS_GNU_MathProg).

If you think you have identified a bug, or want to suggest an enhancement, please log an issue 
[here](https://github.com/OSeMOSYS/OSeMOSYS_GNU_MathProg/issues/new/choose).

### Training materials

- [Understanding the Energy System Lecture Slides](http://www.osemosys.org/understanding-the-energy-system.html) 
by OSeMOSYS Team
- [OSeMOSYS Video Lectures](https://www.youtube.com/watch?v=U9Z4GE_l9mQ&list=PLY5NLA2BTufKU_wDSp-JzP6dhpzwsg1Xx) by OSeMOSYS Team

### OSeMOSYS website

The OSeMOSYS website at [osemosys.org](http://osemosys.org) contains a lot of useful information.

### Newsletter

Sign up for a monthly OSeMOSYS [newsletter](http://www.osemosys.org/news-and-events.html).

## Citing OSemOSYS

If you use OSeMOSYS in your work, please cite the original paper:

    Howells, M., Rogner, H., Strachan, N., Heaps, C., Huntington, H., Kypreos, S., Hughes, A., Silveira, S., DeCarolis, J., Bazilian, M., Roehrl, A., 2011. OSeMOSYS: The Open Source Energy Modeling System: An introduction to its ethos, structure and development. Energy Policy, 39 (10), pp. 5850-5870.

If you use Bibtex, you can use the following:

    @article{HOWELLS20115850,
    title = "OSeMOSYS: The Open Source Energy Modeling System: An introduction to its ethos, structure and development",
    journal = "Energy Policy",
    volume = "39",
    number = "10",
    pages = "5850 - 5870",
    year = "2011",
    note = "Sustainability of biofuels",
    issn = "0301-4215",
    doi = "https://doi.org/10.1016/j.enpol.2011.06.033",
    url = "http://www.sciencedirect.com/science/article/pii/S0301421511004897",
    author = "Mark Howells and Holger Rogner and Neil Strachan and Charles Heaps and Hillard Huntington and Socrates Kypreos and Alison Hughes and Semida Silveira and Joe DeCarolis and Morgan Bazillian and Alexander Roehrl",
    keywords = "Energy systems analysis, Energy modeling, Open source",
    abstract = "This paper discusses the design and development of the Open Source Energy Modeling System (OSeMOSYS). It describes the model’s formulation in terms of a ‘plain English’ description, algebraic formulation, implementation—in terms of its full source code, as well as a detailed description of the model inputs, parameters, and outputs. A key feature of the OSeMOSYS implementation is that it is contained in less than five pages of documented, easily accessible code. Other existing energy system models that do not have this emphasis on compactness and openness makes the barrier to entry by new users much higher, as well as making the addition of innovative new functionality very difficult. The paper begins by describing the rationale for the development of OSeMOSYS and its structure. The current preliminary implementation of the model is then demonstrated for a discrete example. Next, we explain how new development efforts will build on the existing OSeMOSYS codebase. The paper closes with thoughts regarding the organization of the OSeMOSYS community, associated capacity development efforts, and linkages to other open source efforts including adding functionality to the LEAP model."
    }