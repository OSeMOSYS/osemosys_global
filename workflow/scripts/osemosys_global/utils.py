"""Utility functions"""

import pandas as pd
from constants import SET_DTYPES
from typing import Optional

import logging 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

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