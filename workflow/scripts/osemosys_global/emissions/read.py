"""Module for reading in data sources"""

import pandas as pd

from data import _format_ember_emission_data

def import_emission_factors(f: str) -> pd.DataFrame:
    """Imports user defined emission factors.
    
    emission_factors.csv
    """
    return pd.read_csv(f)

def import_iar_base(f: str) -> pd.DataFrame:
    """Imports InputActivityRatio.csv.
    
    InputActivityRatio.csv
    """
    return pd.read_csv(f)

def import_oar_base(f: str) -> pd.DataFrame:
    """Imports OutputActivityRatio.csv.
    
    OutputActivityRatio.csv
    """
    return pd.read_csv(f)

def _read_ember_data(csv_file: str) -> pd.DataFrame:
    """Reads *.csv ember data from https://ember-climate.org/data-catalogue/yearly-electricity-data/

    Data - Yearly Full Release Long Format
    """
    df = pd.read_csv(csv_file)
    df = df.rename(
        columns={"Country code": "COUNTRY", "Year": "YEAR", "Value": "VALUE"}
    )
    df = df[["COUNTRY", "YEAR", "Category", "Subcategory", "Variable", "Unit", "VALUE"]]
    return df[(df.YEAR >= 2015) & (df.Unit != "%")].copy()

def get_ember_emissions(csv_file: str, **kwargs) -> pd.DataFrame:
    df = _read_ember_data(csv_file)
    return _format_ember_emission_data(df)