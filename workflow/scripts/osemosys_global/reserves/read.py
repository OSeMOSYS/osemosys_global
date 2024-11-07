"""Module for reading in data sources"""

import pandas as pd

def import_technologies(f: str) -> pd.DataFrame:
    """Imports TECHNOLOGY.csv as output from the Storage rule."""
    return pd.read_csv(f)

def import_fuels(f: str) -> pd.DataFrame:
    """Imports FUEL.csv as output from the Storage rule."""
    return pd.read_csv(f)