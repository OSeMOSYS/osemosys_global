# Scripts

A collection of useful scripts for running the GNU MathProg implementation of
OSeMOSYS

## cplex_batchRun.bat

This batch script is for running multiple runs in Windows using GNU MathProg
and CPLEX. Place all files together with the batch file, and make sure to change
the naming to your text files and OSeMOSyS version.

To run: just double click the batch file and it will execute.

## kth_hpc_cplex_batchrun.sh

Acknowledgment to Jing Gong at PDC center @ KTH for rewriting the script for Tegner.
This batch script is for running one run for GNU MathProg and CPLEX on High Performance Comupter using SLURM commands
and Linux operating system (such as Tegner @ KTH).
Place all files together with the batch file, and make sure to change the naming
to your text files and OSeMOSyS version.
To run on HPC Tegner @ KTH write: `sbatch kth_hpc_cplex_batchrun.sh`
and it will queue the file and run it.

## convert_cplex_to_cbc.py

This script converts the solution file output of a OSeMOSYS run with the CPLEX
solver into the same format of that produced by CBC.

This removes the zeros, and moves closer to a tidy format that is generally
preferable (sparser and more transferable).

To run use the command `python convert_cplex_to_cbc.py cplex_file output_file`

Use the `-h` or `--help` flag to get the full set of options:

```python
$ python convert_cplex_to_cbc.py --help
usage: convert_cplex_to_cbc.py [-h] [-s START_YEAR] [-e END_YEAR]
                               [--csv | --cbc]
                               cplex_file output_file

Convert OSeMOSYS CPLEX files into different formats

positional arguments:
  cplex_file            The filepath of the OSeMOSYS cplex output file
  output_file           The filepath of the converted file that will be
                        written

optional arguments:
  -h, --help            show this help message and exit
  -s START_YEAR, --start_year START_YEAR
                        Output only the results from this year onwards
  -e END_YEAR, --end_year END_YEAR
                        Output only the results upto and including this year
  --csv                 Output file in comma-separated-values format
  --cbc                 Output file in CBC format, (default option)
```

## preprocess_data.py

This script pre-processes an OSeMOSYS input data file by adding lines that list commodity-technology-mode combinations that data is provided for. Pre-processing a data file before starting a model run significantly reduces the time taken for matrix generation. 

Pre-processing consists of the following steps:
1. Reading the `InputActivityRatio` and `OutputActivityRatio` sections of the data file to identify commodity-technology-mode combinations that data has been explicitly provided for.
2. Adding a set entry for each commodity that lists all technology-mode combinations that are associated with it.  
3. Values from the `InputActivityRatios` and `OutputActivityRatios` sections are added to the sets `MODExTECHNOLOGYperFUELin` and `MODExTECHNOLOGYperFUELout` respectively.
4. Values from the `TechnologyToStorage` and `TechnologyFromStorage` sections are added to the sets `MODExTECHNOLOGYperSTORAGEto` and `MODExTECHNOLOGYperSTORAGEfrom` respectively.
5. All values for technology-mode combinations are added to the sets `MODEperTECHNOLOGY`.

This pre-processing can be run on a terminal window with the following command:
```
python preprocess_data.py <otoole/momani> <input_data_file.txt> <preprocessed_data_file.txt>
``` 

In order to start a model run with a pre-processed data file, the following sets need to be included in the associated OSeMOSYS model file:
```
set MODEperTECHNOLOGY{TECHNOLOGY} within MODE_OF_OPERATION;
set MODExTECHNOLOGYperFUELout{COMMODITY} within MODE_OF_OPERATION cross TECHNOLOGY;
set MODExTECHNOLOGYperFUELin{COMMODITY} within MODE_OF_OPERATION cross TECHNOLOGY;
set MODExTECHNOLOGYperSTORAGEto{STORAGE} within MODE_OF_OPERATION cross TECHNOLOGY;
set MODExTECHNOLOGYperSTORAGEfrom{STORAGE} within MODE_OF_OPERATION cross TECHNOLOGY;
```
## RunHPCNodes.py

This script was written by Taco Niet to run a large number of scenarios on a HPC system managed by the SLURM scheduler.  The script creates the data files for a large number of scenarios using pre-set scenario parameter lists and then submits the list of jobs to full nodes of the HPC system using the parallel command to utilize the full node.  This script works best when there are a large number of scenarios that each take a very short amount of time to run.

First edit the python script for the parameters and scenario combinations you want to investigate.  Then place the 'params_base.dat' file in the directory along with the 'RunHPCNodes.py' script and the OSeMOSYS.mod file.  Log into the login node of the HPC cluster and run the python file.

Note:  The number of scenarios can quickly become quite large with this script so managing the number of runs and how long of a wall time they need is important to consider.  Choose the walltime, number of runs to group and other parameters carefully.
