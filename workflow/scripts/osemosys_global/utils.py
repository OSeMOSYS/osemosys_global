"""Utility Functions"""

import pandas as pd
from typing import Dict, Optional
from pathlib import Path
from osemosys_global.constants import SET_DTYPES

import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

def apply_timeshift(x, timeshift):
    """Applies timeshift to organize dayparts.
    
    Arguments:
        x = Value between 0-24
        timeshift = value offset from UTC (-11 -> +12)"""

    x += timeshift
    if x > 23:
        return x - 24
    elif x < 0:
        return x + 24
    else:
        return x
    
def read_csv(dirpath: str) -> Dict[str,pd.DataFrame]:
    """Reads in CSVs folder
    
    Replace with ReadCSV.read() from otoole v1.0
    """
    data = {}
    files = [Path(x) for x in Path(dirpath).iterdir()]
    for f in files:
        data[f.stem] = pd.read_csv(f)
    return data

def filter_transmission_techs(df: pd.DataFrame, column_name: str = "TECHNOLOGY") -> pd.DataFrame:
    """Filters out only transmission technologies
    
    Arguments: 
        df: pd.DataFrame
            otoole formatted dataframe 
        column_name: str
            Column name to filter on 
            
    Returns: 
        pd.DataFrame
    """
    return df.loc[df[column_name].str.startswith("TRN")].reset_index(drop=True)

def apply_dtypes(df:pd.DataFrame, name: Optional[str]) -> pd.DataFrame:
    """Sets datatypes on dataframe"""
    
    for col in df.columns:
        try:
            df[col] = df[col].astype(SET_DTYPES[col])
        except KeyError:
            if name:
                logging.info(f"Can not set dtype for {name} on {col}")
            else:
                logging.info(f"Can not set dtype on {col}")
    return df
