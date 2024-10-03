"""Function to create sets."""

import pandas as pd

from data import get_years

from constants import SET_DTYPES

def create_sets(x, df, output_dir, custom_node_elements):
    """Creates a formatted otoole set csv 
    
    Arguments: 
        x = Name of set given by otoole as a string
        df = dataframe in extract set information from 
        output_dir = directory to write file to
    
    Returns: 
        None. Writes out csv file  
    
    Example:
        create_sets('TECHNOLOGY', df_oar_final, output_dir)
    """
    set_elements = list(df[x].unique()) + list(df[x].unique()) + list(custom_node_elements)
    set_elements = list(set(set_elements))
    set_elements = [x for x in set_elements if x != 'nan']
    set_elements.sort()
    set_elements_df = pd.DataFrame(set_elements, columns = ['VALUE'])
    return set_elements_df
    
def output_sets(mode_list, start_year, end_year, region_name):                     

    # ## Create set for YEAR, REGION, MODE_OF_OPERATION
    years_df = pd.DataFrame(get_years(start_year, end_year), columns = ['VALUE']).astype(SET_DTYPES["YEAR"])

    mode_list_df = pd.DataFrame(mode_list, columns = ['VALUE']).astype(SET_DTYPES["MODE_OF_OPERATION"])

    regions_df = pd.DataFrame(columns = ['VALUE']).astype(SET_DTYPES["REGION"])
    regions_df.loc[0] = region_name
    
    return years_df, mode_list_df, regions_df