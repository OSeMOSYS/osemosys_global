"""Custom Nodes logic"""

import pandas as pd
import itertools


def _get_custom_demand_expected(
    nodes: list[str], start_year: int, end_year: int
) -> pd.DataFrame:
    """Gets formatted expected custom data"""

    years = range(start_year, end_year + 1)

    df = pd.DataFrame(
        list(itertools.product(nodes, years)), columns=["CUSTOM_NODE", "YEAR"]
    )
    df["REGION"] = "GLOBAL"
    df["FUEL"] = "ELC" + df["CUSTOM_NODE"] + "02"

    return df


def import_custom_demand_data(csv: str) -> pd.DataFrame:
    """Gets all custom demand data"""
    return pd.read_csv(csv)


def get_custom_demand_data(
    all_custom: pd.DataFrame, start_year: int, end_year: int
) -> pd.DataFrame:
    """Gets merged custom demand data"""
    
    nodes = all_custom['CUSTOM_NODE'].unique()
    expected = _get_custom_demand_expected(nodes, start_year, end_year)

    df = pd.merge(expected, all_custom, how="left", on=["CUSTOM_NODE", "YEAR"])
    df = df[["REGION", "FUEL", "YEAR", "VALUE"]]

    return df


def merge_default_custom_data(
    default: pd.DataFrame, custom: pd.DataFrame
) -> pd.DataFrame:
    assert default.columns.equals(custom.columns)
    df = pd.concat([default, custom], ignore_index=True)
    df["VALUE"] = df["VALUE"].round(2)
    df = df.drop_duplicates(keep="last", subset=["REGION", "FUEL", "YEAR"])
    return df
