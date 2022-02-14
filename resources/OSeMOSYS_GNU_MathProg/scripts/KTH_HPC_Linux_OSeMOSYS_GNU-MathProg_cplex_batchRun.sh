######### High Performance Computer (HPC) specific commands (SLURM) ########

#!/bin/bash 

# Set the allocation to be charged for this job
#SBATCH -A 

# The name of the script is myjob
#SBATCH -J myjob

# 10 hour wall-clock time will be given to this job
#SBATCH -t 10:00:00

# Number of nodes
#SBATCH --nodes=1
# Number of MPI processes per node
#SBATCH --ntasks-per-node=24

#SBATCH --mail-type=ALL

#SBATCH -e error_file.e
#SBATCH -o output_file.o

############### Here finished the SLURM commands #############

# Load the glpk mode which includes glpsol, please check version this command is specific for Tegner
# Set the path to glpk (need to change) if not on Tegner
module add glpk/4.65

# Set the path to cplex here

# this command starts glpsol (GNU Mathprog) 
# Only generate the .lp file using flag --check

glpsol -m OSeMOSYS.txt -d data.txt --wlp matrix.lp --check

# break mean make a new empty file for mycplexcommands
rm -f mycplexcommands
touch mycplexcommands

# echo writes each line to mycplexcommands that I want to execute in CPLEX

echo "read matrix.lp" 		> mycplexcommands
echo "optimize"             >> mycplexcommands
echo "write"                >> mycplexcommands
echo "matrix.sol"    		>> mycplexcommands
echo "quit"                 >> mycplexcommands


# xecutes the cplex script written above
# Should set the path to CPLEX
cplex < mycplexcommands

# the sol file is input to transform python script
#
python transform_31072013.py matrix.sol matrix.txt


