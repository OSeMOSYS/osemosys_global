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

###
# public functions
###


def get_eia_capacity(json_file: str) -> pd.DataFrame:
    df = _read_eia_data(json_file)
    return _format_eia_capacity_data(df)


def get_eia_generation(json_file: str) -> pd.DataFrame:
    df = _read_eia_data(json_file)
    return _format_eia_generation_data(df)


def format_og_generation(prod_tech_annual: pd.DataFrame) -> pd.DataFrame:
    """Formats ProductionByTechnologyAnnual data for eia comparison"""

    name_mapper = {
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

    df = prod_tech_annual.copy()

    if len(df.columns) == 1:
        df = df.reset_index()

    df = df[(df.TECHNOLOGY.str.startswith("PWR")) & (df.YEAR < 2023)]
    df["COUNTRY"] = df.TECHNOLOGY.str[6:9]
    df["CODE"] = df.TECHNOLOGY.str[3:6]
    df["CODE"] = df.CODE.map(name_mapper)
    df = df.dropna(subset="CODE")
    df["TECHNOLOGY"] = df.CODE + df.COUNTRY
    df = df.drop(columns=["FUEL", "COUNTRY", "CODE"])
    return df.groupby(["REGION", "TECHNOLOGY", "YEAR"]).sum()


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
    df["year"] = df.data.map(lambda x: datetime.fromtimestamp(x["date"] / 1000).year)
    df["VALUE"] = df.data.map(lambda x: x["value"])
    df["VALUE"] = df.VALUE.fillna(0)
    return df.drop(
        columns=["series_id", "frequency", "productid", "activityid", "data", "unit"]
    )


def _format_eia_capacity_data(eia: pd.DataFrame) -> pd.DataFrame:
    """Formats data into otoole compatiable data structure"""

    df = eia.copy()

    df["name"] = df.name.map(
        lambda x: x.split(" electricity installed capacity")[0]
    ).map(CAPACITY_MAPPER)
    df["name"] = df.name + df.iso
    df = df.drop(columns=["iso"])
    df = df.groupby(["name", "year"], as_index=False).sum()
    df = df.rename(columns={"name": "TECHNOLOGY", "year": "YEAR"})
    df["REGION"] = "GLOBAL"
    return df.set_index(["REGION", "TECHNOLOGY", "YEAR"])


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
    df = df.groupby(["name", "year"], as_index=False).sum()
    df = df.rename(columns={"name": "TECHNOLOGY", "year": "YEAR"})
    df["REGION"] = "GLOBAL"
    # billion kWh -> PJ
    # 1B kWh = 1 TWh * (1PWh / 1000TWh) * (3600sec / hr) = 1 PWs = 1 PJ
    df["VALUE"] = df.VALUE.mul(3.6)
    return df.set_index(["REGION", "TECHNOLOGY", "YEAR"])
