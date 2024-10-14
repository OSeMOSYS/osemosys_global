"""Data handeling for Climate Watch validation

https://www.climatewatchdata.org/ghg-emissions
"""

import pandas as pd

###
# public functions
###


def get_cw_emissions(csv_file: str, **kwargs) -> pd.DataFrame:
    """Gets Climate Watch emissions data"""
    df = _read_cw_data(csv_file)
    return _format_cw_data(df)


###
# private functions
###


def _read_cw_data(csv_file: str) -> pd.DataFrame:
    """Reads climate watch data

    https://www.climatewatchdata.org/ghg-emissions?end_year=2021&gases=all-ghg&sectors=electricity-heat&start_year=1990
    """
    return pd.read_csv(csv_file, skipfooter=2, engine="python")


def _format_cw_data(cw: pd.DataFrame) -> pd.DataFrame:
    df = cw.copy()
    df = df.drop(columns=["Country/Region", "unit"]).rename(columns={"iso": "EMISSION"})
    df = df.melt(id_vars=["EMISSION"], var_name="YEAR", value_name="VALUE")
    df = df.fillna(0).replace("false", 0)
    df["YEAR"] = df.YEAR.astype(int)
    df["VALUE"] = df.VALUE.astype(float)
    df["REGION"] = "GLOBAL"
    return df.set_index(["REGION", "EMISSION", "YEAR"])
