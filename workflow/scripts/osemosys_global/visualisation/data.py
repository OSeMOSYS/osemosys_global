"""Getter functions for data to plot"""

from .utils import powerplant_filter, transform_ts
import pandas as pd
from typing import Dict
import logging 
logger = logging.getLogger(__name__)

def get_total_capacity_data(data: Dict[str,pd.DataFrame], country:str =None):
    """ Gets data for plotting total capacity 
    
    Arguments: 
        data: Dict[str,pd.DataFrame]
            Internal datastore 
        country: str 
            If a country provided, plot at a country level, else plot at a 
            system level
    
    Returns
        df: pd.DataFrame 
            Dataframe formatted for total capacity data 
    """
    df = data["TotalCapacityAnnual"]
    df = powerplant_filter(df, country)
    df.VALUE = df.VALUE.astype('float64')
    df = df.groupby(['LABEL', 'YEAR'],
                    as_index=False)['VALUE'].sum()
    return df

def get_generation_annual_data(data: Dict[str,pd.DataFrame], country:str =None):
    """Gets data for plotting annual generation
    
    Arguments: 
        data: Dict[str,pd.DataFrame]
            Internal datastore 
        country: str 
            If a country provided, plot at a country level, else plot at a 
            system level
    
    Returns 
        df: pd.DataFrame 
            Dataframe formatted for total annual generation
    """
    df = data["ProductionByTechnologyAnnual"]
    df = powerplant_filter(df, country)
    df.VALUE = df.VALUE.astype('float64')
    df = df.groupby(['LABEL', 'YEAR'],
                    as_index=False)['VALUE'].sum()
    return df

def get_generation_ts_data(input_data: Dict[str,pd.DataFrame], result_data: Dict[str,pd.DataFrame], country:str =None):
    """Gets data for plotting generation by time slice
    
    Arguments: 
        input_data: Dict[str,pd.DataFrame]
            Internal datastore for input data
        result_data: Dict[str,pd.DataFrame]
            Internal datastore for results
        country: str 
            If a country provided, plot at a country level, else plot at a 
            system level
    
    Returns 
        df: pd.DataFrame 
            Dataframe formatted for total annual generation"""
            
    df = result_data["ProductionByTechnology"]
    df = powerplant_filter(df, country)
    df.VALUE = df.VALUE.astype('float64')
    df = transform_ts(input_data, df)
    return df