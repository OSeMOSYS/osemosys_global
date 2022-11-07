"""Extracts results from the list of input files

Notes
-----
Note that this script is directly implemented into the Snakemake
workflow and pulls the arguments for the ``main()`` function directly from
``snakemake.input`` and ``snakemake.output[0]`` attributes on the snakemake
object passed into this module at run time.
"""

import pandas as pd
from typing import List, Tuple, Dict
from pathlib import Path
from utils import get_model_run_scenario_from_filepath, parse_yaml
import sys
from utils import write_results, read_results
import itertools

from logging import getLogger

logger = getLogger(__name__)

def get_indices(parameters : pd.DataFrame, filename : str) -> Dict :
    indices = parameters.set_index('filename').loc[filename].dropna().drop('resultfile')
    indices = indices.to_dict()
    return {x:indices[x].split(',') for x in indices}

def main(input_files: List, output_file: str, indices: Tuple, config: Dict):
    """Iterate over list of CSV files, extract defined results, write to output file.

    Parameters
    ----------
    input_files : List
        List of CSVs to iterate over
    output_file : str
        Name of output file
    indices : Dict
        Indices to extract value over 
        ex. {'REGION':['SIMPLICITY], 'TECHNOLOGY':['GAS_EXTRACTION','HYD1']}
    config : Dict
        Datapackage holding type information 
    """
    aggregated_results = []
    for filename in input_files:

        bits = get_model_run_scenario_from_filepath(filename)
        df_index = config[bits['param']]['indices']
        column_dtypes = {
            x:config[x]['dtype'] for x in df_index
        }
        df = read_results(filename)
        df = df.astype(column_dtypes).set_index(df_index)
        ################################################################
        # this method of slicing and then appending is super enefficient... 
        # will revist to speed this up 
        result_dfs = []
        parameters = tuple(itertools.product(*indices.values()))
        indices_expanded = tuple([tuple(indices.keys())] * len(parameters))
        # print(parameters)
        # print(indices_expanded)
        for index, param in zip(indices_expanded, parameters):
            try:
                results = df.xs(param, level=index, drop_level=False)
                result_dfs.append(results)
            except KeyError as ex:
                raise ex
        results = pd.concat(result_dfs)
        results = results.reset_index(level='YEAR')
        ################################################################
        # results['SCENARIO'] = bits['scenario']
        results['MODELRUN'] = bits['model_run']
        # results = results.reset_index(
        #     ).set_index(['SCENARIO', 'MODELRUN'] + df_index)
        results = results.reset_index(
            ).set_index(['MODELRUN'] + df_index)
        aggregated_results.append(results)

    results = pd.concat(aggregated_results)
    write_results(results, output_file, True)

if __name__ == '__main__':

    if "snakemake" in globals():
        input_files = snakemake.input['csvs']
        yaml_config = snakemake.input['config']
        output_file_path = snakemake.output[0]
        indices = snakemake.params['parameter']
    else:
        if len(sys.argv) != 5:
            raise ValueError(
                "Usage: python extract_results.py <input_csvs> <user_config> <output_file> <result_parameters>"
            )
        input_files = sys.argv[1]
        if not isinstance(input_files, list):
            input_files = [input_files]
        yaml_config = sys.argv[2]
        output_file_path = sys.argv[3]
        parameters = pd.read_csv(sys.argv[4])
        output_file = Path(sys.argv[3]).stem
        indices = get_indices(parameters, output_file)

    if 'YEAR' in indices:
        indices['YEAR'] = [float(x) for x in indices['YEAR']]

    user_config = parse_yaml(yaml_config)
    main(input_files, output_file_path, indices, user_config)