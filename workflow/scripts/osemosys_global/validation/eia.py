"""Data handeling for EIA validation"""

import pandas as pd
from datetime import datetime

###
# constants for data allignment
###

CAPACITY_MAPPER = {
    "Hydroelectric pumped storage": "HPS",
    "Tide and wave": "WAV",
    "Nuclear": "URN",
    "Non-hydro renewable": "NHR",
    "Hydroelectricity installed capacity": "HYD",
    "Fossil fuels": "FFS",
    "Renewable": "RNW",
    "Biomass and waste": "BIO",
    "Wind": "WND",
    "Solar": "SPV",
    "Geothermal": "GEO",
    "Electricity installed capacity": "ELC",
    "Renewable": "RNW",
}

GENERATION_MAPPER = {
    "Hydroelectric pumped storage": "HPS",
    "Tide and wave": "WAV",
    "Nuclear": "URN",
    "Non-hydro renewable": "NHR",
    "Hydroelectricity net generation": "HYD",
    "Fossil fuels": "FFS",
    "Renewable": "RNW",
    "Biomass and waste": "BIO",
    "Wind": "WND",
    "Solar": "SPV",
    "Geothermal": "GEO",
    "Electricity net generation": "ELC",
    "Renewable": "RNW",
    "Petroleum fossil fuel": "OIL",
    "Coal fossil fuel": "COA",
    "Natural gas fossil fuel": "GAS",
    "Other gases fossil fuel": "OTH",
}

OG_GEN_NAME_MAPPER = {
    "BIO": "BIO",
    "CCG": "GAS",
    "COA": "COA",
    "CSP": "SPV",
    "HYD": "HYD",
    "OCG": "GAS",
    "OIL": "OIL",
    "SPV": "SPV",
    "TRN": None,
    "URN": "URN",
    "WON": "WND",
    "WOF": "WND",
    "WAV": "WAV",
}

OG_CAP_NAME_MAPPER = {
    "BIO": "BIO",
    "CCG": "FFS",
    "COA": "FFS",
    "CSP": "SPV",
    "HYD": "HYD",
    "OCG": "FFS",
    "OIL": "FFS",
    "SPV": "SPV",
    "TRN": None,
    "URN": "URN",
    "WON": "WND",
    "WOF": "WND",
    "WAV": "WAV",
}

###
# public functions
###


def get_eia_capacity(json_file: str, **kwargs) -> pd.DataFrame:
    df = _read_eia_data(json_file)
    return _format_eia_capacity_data(df)


def get_eia_generation(json_file: str, **kwargs) -> pd.DataFrame:
    df = _read_eia_data(json_file)
    return _format_eia_generation_data(df)


###
# private functions
###


def _read_eia_data(json_file: str) -> pd.DataFrame:
    """Reads *.json EIA data from https://www.eia.gov/international/data/world

    Data -> Electricity -> Electricity Capacity / Electricity Generation
    """
    df = pd.read_json(json_file)
    df["name"] = df.name.map(lambda x: x.split(", ")[0])
    df = df.explode(column="data")
    # not sure why, but the 'datetime.fromtimestamp(x["date"] / 1000).year' call gives the
    # next year rather than the correct one. ie. If I call the year 2020, 2021 values are
    # returned. Thats why the extra '+1' at the end of the lambda
    df["year"] = df.data.map(
        lambda x: datetime.fromtimestamp(x["date"] / 1000).year + 1
    )
    df["VALUE"] = df.data.map(lambda x: x["value"])
    df["VALUE"] = df.VALUE.fillna(0)
    return df.drop(
        columns=["series_id", "frequency", "productid", "activityid", "data", "unit"]
    )


def _format_eia_capacity_data(eia: pd.DataFrame) -> pd.DataFrame:
    """Formats data into otoole compatiable data structure

    Note, no unit conversion as capacity is already given in GW
    """

    df = eia.copy()

    df["name"] = df.name.map(
        lambda x: x.split(" electricity installed capacity")[0]
    ).map(CAPACITY_MAPPER)
    df["name"] = df.name + df.iso
    df["VALUE"] = (
        df.VALUE.replace("NA", "0")
        .replace("--", "0")
        .replace("ie", "0")
        .replace("(s)", "0")
        .fillna("0")
        .astype(float)
    )
    df = df.rename(columns={"name": "TECHNOLOGY", "year": "YEAR"})
    df["REGION"] = "GLOBAL"
    df = df[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    return df.groupby(["REGION", "TECHNOLOGY", "YEAR"]).sum()


def _format_eia_generation_data(eia: pd.DataFrame) -> pd.DataFrame:
    """Formats data into otoole compatiable data structure"""

    df = eia.copy()

    df["name"] = df.name.map(lambda x: x.split(" electricity net generation")[0]).map(
        GENERATION_MAPPER
    )
    df["name"] = df.name + df.iso
    df = df.drop(columns=["iso"])
    df["VALUE"] = (
        df.VALUE.replace(r"^se\|", "", regex=True)
        .replace("ie", 0)
        .replace("(s)", 0)
        .replace("--", 0)
        .replace("NA", 0)
        .astype(float)
    )
    df = df.rename(columns={"name": "TECHNOLOGY", "year": "YEAR"})
    df["REGION"] = "GLOBAL"
    # billion kWh -> PJ
    # 1B kWh = 1 TWh * (1PWh / 1000TWh) * (3600sec / hr) = 1 PWs = 1 PJ
    df["VALUE"] = df.VALUE.mul(3.6)
    df = df[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    return df.groupby(["REGION", "TECHNOLOGY", "YEAR"]).sum()
