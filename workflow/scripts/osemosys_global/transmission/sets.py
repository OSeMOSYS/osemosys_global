"""Function to create sets."""

import pandas as pd

def get_unique_technologies(ar: pd.DataFrame) -> list[str]:
    """Gets unique technologies from activity ratio dataframe"""
    return ar.TECHNOLOGY.unique().tolist()

def get_unique_fuels(ar: pd.DataFrame) -> list[str]:
    """Gets unique fuels from activity ratio dataframe"""
    return ar.FUEL.unique().tolist()

def create_set_from_iterators(*args: list[str] | set[str]) -> pd.DataFrame:
    """Creates a set formatted dataframe from an arbitrary number of arguments"""
    
    data = []
    
    for arg in args:
        data += arg
    
    data = list(set(data))
    
    return pd.Series(data, name="VALUE").to_frame().sort_values(['VALUE'])