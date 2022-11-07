"""Generates an unscaled sample from a list of parameters

Arguments
---------
path_to_parameters : str
    File containing the parameters to generate a sample for

sample_file : str
    File path to save sample to

replicates : int
    The number of model runs to generate

Usage
-----
To run the script on the command line, type::

    python create_sample.py path/to/parameters.csv path/to/save.txt 10

The ``parameters.csv`` CSV file should be formatted as follows::

    name,group,indexes,min_value,max_value,dist,interpolation_index,action
    CapitalCost,pvcapex,"GLOBAL,GCPSOUT0N",500,1900,unif,YEAR,interpolate
    DiscountRate,discountrate,"GLOBAL,GCIELEX0N",0.05,0.20,unif,None,fixed

"""
from SALib.sample import morris
import os
import numpy as np
import csv
from typing import List
import sys
import utils

from logging import getLogger

logger = getLogger(__name__)

def main(parameters: dict, sample_file: str, replicates: int):

    problem = utils.create_salib_problem(parameters)

    sample = morris.sample(problem, N=50, optimal_trajectories=replicates,
                           local_optimization=True, seed=42)

    np.savetxt(sample_file, sample, delimiter=',')


if __name__ == "__main__":

    parameters_file = sys.argv[1]
    sample_file = sys.argv[2]
    replicates = int(sys.argv[3])
    with open(parameters_file, 'r') as csv_file:
        reader = list(csv.DictReader(csv_file))

    main(reader, sample_file, replicates)
