import os
import pandas as pd
from typing import List, Dict
import yaml
from pathlib import Path

from logging import getLogger

logger = getLogger(__name__)

def get_model_run_scenario_from_input_filepath(filename: str):
    """Parses filepath to extract useful bits

    "results/{scenario}/SA/modelruns/model_{modelrun}/data/{input_file}.csv"
    """

    path = Path(filename)
    filepath = path.parent
    param = path.with_suffix('').name
    scenario = path.parts[1]
    model_run = path.parts[4]
    return {'model_run': model_run, 'scenario': scenario,
            'param': param, 'filepath': filepath}

def get_model_run_scenario_from_filepath(filepath: str) -> Dict:
    """Parses filepath to extract useful bits.

    Input filepath is expected in the form of:
        "results/{scenario}/SA/modelruns/model_{modelrun}/results/{input_file}.csv"

    Parameters
    ----------
    file : str
        file path from root directory 

    Returns
    -------
    Dict
        With model_run, scenario, OSeMOSYS parameter, and filepath directory
    """
    f = Path(filepath)
    parts = f.parts
    return {'model_run': parts[4].split('_')[1], 'scenario': parts[1],
            'param': f.stem, 'filepath': str(f.parent)}

def read_results(input_filepath: str) -> pd.DataFrame:
    extension = os.path.splitext(input_filepath)[1]
    if extension == '.parquet':
        df = pd.read_parquet(input_filepath)
    elif extension == '.csv':
        df = pd.read_csv(input_filepath)
    elif extension == '.feather':
        df = pd.read_feather(input_filepath)
    return df


def write_results(df: pd.DataFrame, output_filepath: str, index=None) -> None:
    """Write out aggregated results to disk by scenario 

    Arguments
    ---------
    df: pd.DataFrame
        Dataframe to write out
    output_filepath: str
        Path to the output file
    index=None
        Whether to write out the index or not
    """
    extension = os.path.splitext(output_filepath)[1]
    if extension == '.parquet':
        df.to_parquet(output_filepath, index=index)
    elif extension == '.csv':
        df.to_csv(output_filepath, index=index)
    elif extension == '.feather':
        if index:
            df = df.reset_index()
        df.to_feather(output_filepath)

def create_salib_problem(parameters: List) -> dict:
    """Creates SALib problem from scenario configuration.
    
    Arguments
    ---------
    parameters: List
        List of dictionaries describing problem. Each dictionary must have
        'name', 'indexes', 'group' keys

    Returns
    -------
    problem: dict
        SALib formatted problem dictionary

    Raises
    ------
    ValueError
        If only one variable is givin, OR 
        If only one group is given
    """

    problem = {}
    problem['num_vars'] = len(parameters)
    if problem['num_vars'] <= 1:
        raise ValueError(
            f"Must define at least two variables in problem. User defined "
            f"{problem['num_vars']} variable(s).")

    names = []
    bounds = []
    groups = []
    for parameter in parameters:
        names.append(parameter['name'] + ";" + parameter['indexes'])
        groups.append(parameter['group'])
        min_value = 0
        max_value = 1
        bounds.append([min_value, max_value])

    problem['names'] = names
    problem['bounds'] = bounds
    problem['groups'] = groups
    num_groups = len(set(groups))
    if num_groups <= 1:
        raise ValueError(
            f"Must define at least two groups in problem. User defined "
            f"{num_groups} group(s).")

    return problem

def parse_yaml(path : str) -> Dict:
    """Parses a YAML file to a dictionary

    Parameters
    ----------
    path : str
        input path the yaml file 

    Returns
    -------
    parsed_yaml : Dict
        parsed YAML file

    Raises
    ------
    YAMLError
        If the yaml file can't be loaded 
    """
    with open(path, 'r') as stream:
        try:
            parsed_yaml = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
        return parsed_yaml