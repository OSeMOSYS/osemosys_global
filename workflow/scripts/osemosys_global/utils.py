"""Utility Functions"""

import pandas as pd
from typing import Optional
from constants import SET_DTYPES

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
