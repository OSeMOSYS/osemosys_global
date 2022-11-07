"""Creates model runs from a sample file and master model datapackage

Arguments
---------
<input_filepath>
    Path to the master model datapackage
<output_filepath>
    Path to the new model file directory
<sample_filepath>
    Path the sample file

The expected format of the input sample file is a CSV file with the following structure::

    name,indexes,value_base_year,value_end_year,action,interpolation_index
    CapitalCost,"GLOBAL,GCPSOUT0N",1024.3561663863075,2949.23939,interpolate,YEAR

It is very similar to the overall ``parameter.csv`` configuration file, except holds a sample value
rather than a range

To run this script on the command line, use the following::

    python create_modelrun.py <input_filepath> <output_filepath> <sample_filepath>

"""
import sys

from typing import Dict, List, Union, Tuple, Any, Optional

import csv
import pandas as pd
import numpy as np
from otoole.read_strategies import ReadDatapackage
from otoole.write_strategies import WriteCsv
from otoole.utils import _read_file
import os

from logging import getLogger

logger = getLogger(__name__)
def process_data(
                 start_year_value: float,
                 end_year_value: float,
                 first_year: int,
                 last_year: int
                 ) -> np.array:
    """Interpolate data between min and max years

    Arguments
    ---------
    start_year_value: float
        Value of the parameter in the start year
    end_year_value: float
        Value of the parameter in the end year
    first_year: int
        First year of the range to interpolate
    last_year: int
        Last year of the range to interpolate

    Returns
    -------
    values: np.array
        Linearlly interpolated values between the start and end year. 
    """
    values = np.interp([range(int(first_year), int(last_year) + 1)],
                       np.array([int(first_year), int(last_year)]),
                       np.array([float(start_year_value), float(end_year_value)])).T
    return values


def get_types_from_tuple(index: list, param: str, config: Dict) -> Tuple:
    """Checks datatype of a parameters indices.

    Parameters
    ----------
    index : list
        List of indices for the parameters.
    param : str
        Name of parameter
    config : Dict
        Input datapackage

    Returns
    -------
    typed_index : Tuple
        Parameter with properly typed value
    """
    depth = len(index)
    names = config[param]['indices'][0:depth]
    typed_index = []
    dtypes = config[param]['index_dtypes']
    for name, element in zip(names, index):
        this_type = dtypes[name]
        if this_type == 'str':
            typed_index.append(str(element))
        elif this_type == 'float':
            typed_index.append(float(element))
        elif this_type == 'int':
            typed_index.append(int(element))

    return typed_index

def apply_interpolation_action(
    action : str,
    inter_index: str,
    start_year_value : float = None, 
    end_year_value : Optional[float] = None,
    first_year : Optional[int] = None,
    end_year : Optional[int] = None
) -> Any:
    """Applies the interploation method. 

    Parameters
    ----------
    action : str
        'Fixed' or 'interpolate'
    inter_index: str
        'YEAR' or None
    start_year_value : float = None
        Value of the parameter in the start year
    end_year_value : Optional[float] = None,
        Value of the parameter in the end year. Only needed if inter_index is 'YEAR' and action is 'interpolate'
    first_year : Optional[float] = None,
        First year of the range to interpolate. Only needed if inter_index is 'YEAR'
    end_year : Optional[float] = None
        Last year of the range to interpolate. Only needed if inter_index is 'YEAR'

    Returns
    -------
    model_params : Any
        An array with interpolated values over each year, or a single value in a list

    Raises
    ------
    ValueError
        If `action` is not either `fixed` or `interpolate`
    ValueError
        If `inter_index` is not `None` or `YEAR`
    """
    if action == 'interpolate':
        new_values = process_data(start_year_value, end_year_value, first_year, end_year)
    elif action == 'fixed':
        if not inter_index: # None
            new_values = [start_year_value]
        elif inter_index == 'YEAR':
            new_values = np.full((end_year + 1 - first_year, 1), start_year_value)
        else:
            logger.error(f"{inter_index} is not a valid index for the {action} interpolation action")
            raise ValueError
    else:
        logger.error(f"{action} not a valid interpolation action. Choose between 'fixed' and 'interpolate'")
        raise ValueError
    return new_values

def apply_yearly_interploted_values(
    name : str,
    model_params : Dict[str, pd.DataFrame],
    index : Tuple,
    interpolated_values : np.array, 
    first_year : int,
    end_year : int
) -> Dict[str, pd.DataFrame]:
    """Applies interpolated values over years. 

    Parameters
    ----------
    name : str
        name of OSeMOSYS parameter
    model_params : Dict[str, pd.DataFrame]
        Original input OSeMOSYS parameters 
    index : Tuple
        Value of the parameter in the start year
    interpolated_values : np.array
        Array that will replace existing values
    first_year : int
        First year of the range to interpolate
    end_year : int
        Last year of the range to interpolate

    Returns
    -------
    model_params : Dict[str, pd.DataFrame]
        Updated input OSeMOSYS parameters 

    Raises
    ------
    ValueError
        If index not found in original dataframe 
    """
    print(model_params[name])
    try:
        model_params[name].loc[tuple(index + [first_year]):tuple(index + [end_year])] = interpolated_values
    except ValueError as ex:
        logger.error(f"Error raised in parameter {name} by index {[tuple(index + [first_year]), tuple(index + [end_year])]}")
        raise ValueError(ex)
    return model_params

def apply_interploted_values(
    name : str,
    model_params : Dict[str, pd.DataFrame],
    index : Tuple,
    interpolated_values : list, 
) -> Dict[str, pd.DataFrame]:
    """Applies interpolated values that ARE NOT index over a set. 

    Parameters
    ----------
    name : str
        name of OSeMOSYS parameter
    model_params : Dict[str, pd.DataFrame]
        Original input OSeMOSYS parameters 
    index : Tuple
        Value of the parameter in the start year
    interpolated_values : list
        Single item list of value to replce existing value

    Returns
    -------
    model_params : Dict[str, pd.DataFrame]
        Updated input OSeMOSYS parameters 

    Raises
    ------
    ValueError
        If index not found in original dataframe 
    """
    try:
        model_params[name].loc[tuple(index), 'VALUE'] = interpolated_values[0]
    except ValueError as ex:
        logger.error(f"Error raised in parameter {name} by index {index}")
        raise ValueError(ex)
    return model_params

def modify_parameters(
        model_params: Dict[str, pd.DataFrame],
        parameters: List[Dict[str, Union[str, int, float]]],
        config: Dict) ->  Dict[str, pd.DataFrame]:
    """Modifies model parameters based on Morris Sample. 

    The action and interpolation method are read in from the parameters argument
    and applied to the model_params based on the interpolation index. 

    If the interpolation index is set to `YEAR` AND the interpolation method is 
    set to `fixed`, the base year value is set for all years. 

    If the interpolation index is set to `YEAR` AND the interpolation method is 
    set to `interpolate`, the values between the base and end years are linearlly
    interpolated. 

    If the interpolation index is set to `None`, then no interpolation takes place. 

    Parameters
    ----------
    model_params: Dict[str, pd.DataFrame]
        Input OSeMOSYS parameters 
    parameters : List[Dict[str, Union[str, int, float]]]
        Flattened input parameters for the individual model run
    config : Dict
        Input datapackage

    Returns
    -------
    model_params : Tuple
        Parameter with properly typed value
    """

    first_year = model_params['YEAR'].min().values[0]
    end_year = model_params['YEAR'].max().values[0]

    for parameter in parameters:

        # retrieve input data
        name = parameter['name']
        untyped_index = parameter['indexes'].split(",")
        index = get_types_from_tuple(untyped_index, name, config)
        start_year_value = float(parameter['value_base_year'])
        end_year_value = float(parameter['value_end_year'])
        action = parameter['action']
        inter_index = parameter['interpolation_index']

        # retrieve interpolated values
        new_values = apply_interpolation_action(action, inter_index, start_year_value,
            end_year_value, first_year, end_year)

        # apply interpolated values
        logger.info("Updating values for {} in {}".format(index, name))
        if inter_index == 'YEAR':
            model_params = apply_yearly_interploted_values(name, model_params, index,
                new_values, first_year, end_year)
        else:
            model_params = apply_interploted_values(name, model_params, index, new_values)

    return model_params

def main(
    input_filepath, 
    output_filepath, 
    parameters: List[Dict[str, Union[str, int, float]]],
    user_config):

    model_params, default_values = ReadDatapackage(user_config=user_config).read(input_filepath)

    logger.info("Reading datapackage {}".format(input_filepath))
    for name, parameter in model_params.items():
        parameter = parameter.sort_index()
        model_params[name] = parameter
    model_params = modify_parameters(model_params, parameters, user_config)
    # WriteDatapackage(user_config=user_config).write(model_params, output_filepath, default_values)
    WriteCsv(user_config=user_config).write(model_params, output_filepath, default_values)


if __name__ == "__main__":

    if len(sys.argv) != 5:
        print("Usage: python create_modelrun.py <input_filepath> <output_filepath> <sample_filepath> <user_config>")
    else:
        with open(sys.argv[3], 'r') as csv_file:
            sample = list(csv.DictReader(csv_file))

        _, ending = os.path.splitext(sys.argv[4])
        with open(sys.argv[4], "r") as f:
            user_config = _read_file(f, ending)

        main(sys.argv[1], sys.argv[2], sample, user_config)
