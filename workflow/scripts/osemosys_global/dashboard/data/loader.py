import pandas as pd
from otoole import ReadCsv
import geopandas as gpd
from typing import Dict

def load_osemosys_sets(input_data: Dict[str,pd.DataFrame]) -> pd.DataFrame:
    """Loads all set data into single dataframe"""

    df = pd.DataFrame()
    sets = [
        "DAILYTIMEBRACKET",
        "DAYTYPE",
        "EMISSION",
        "FUEL",
        "MODE_OF_OPERATION",
        "REGION",
        "SEASON",
        "STORAGE",
        "TECHNOLOGY",
        "TIMESLICE",
        "YEAR",
    ]
    for s in sets:
        df[s] = input_data[s]["VALUE"]
    return df