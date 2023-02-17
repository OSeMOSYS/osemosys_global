import pandas as pd
from otoole import ReadCsv
import geopandas as gpd
from ...visualisation.utils import parse_node_data_codes, parse_line_data_codes

def load_osemosys_data(csv: str) -> dict[str:pd.DataFrame]:
    """Loads all osemosys input data"""
    input_data, _ = ReadCsv().read(filepath=csv)
    input_data["Demand"] = input_data["SpecifiedAnnualDemand"].mul(input_data["SpecifiedDemandProfile"], fill_value=0.0)
    return input_data

def load_osemosys_sets(csv: str) -> pd.DataFrame:
    """Loads all set data into single dataframe"""
    input_data = load_osemosys_data(csv)
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

def load_node_data(csv: str) -> gpd.GeoDataFrame:
    input_data = pd.read_csv(csv)
    df = parse_node_data_codes(input_data)
    return gpd.GeoDataFrame(
        df, geometry=gpd.points_from_xy(df["Longitude"], df["Latitude"])
    )
    
def load_line_data(csv: str) -> pd.DataFrame:
    input_data = pd.read_csv(csv)
    return parse_line_data_codes(input_data)