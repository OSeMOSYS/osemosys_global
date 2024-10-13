"""Processes data from ember database

https://ember-climate.org/data-catalogue/yearly-electricity-data/
"""

import pandas as pd

###
# constants for data allignment
###

TECHNOLOGY_MAPPER = {
    "Bioenergy": "BIO",
    "Coal": "COA",
    "Gas": "GAS",
    "Hydro": "HYD",
    "Nuclear": "URN",
    "Other Fossil": "OTH",
    "Other Renewables": "OTH",
    "Solar": "SPV",
    "Wind": "WND",
}

OG_NAME_MAPPER = {
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

###
# public functions
###


def get_ember_capacity(csv_file: str, **kwargs) -> pd.DataFrame:
    df = _read_ember_data(csv_file)
    return _format_ember_capacity_data(df)


def get_ember_generation(csv_file: str, **kwargs) -> pd.DataFrame:
    df = _read_ember_data(csv_file)
    return _format_ember_generation_data(df)


def format_og_data(og: pd.DataFrame) -> pd.DataFrame:
    """Formats OG results for ember comparison

    Works on:
    - ProductionByTechnologyAnnual
    - TotalCapacityAnnual
    """

    df = og.copy()

    if len(df.columns) == 1:
        df = df.reset_index()

    df = df[(df.TECHNOLOGY.str.startswith("PWR")) & (df.YEAR < 2023)]
    df["COUNTRY"] = df.TECHNOLOGY.str[6:9]
    df["CODE"] = df.TECHNOLOGY.str[3:6]
    df["CODE"] = df.CODE.map(OG_NAME_MAPPER)
    df = df.dropna(subset="CODE")
    df["TECHNOLOGY"] = df.CODE + df.COUNTRY
    df = df[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    return df.groupby(["REGION", "TECHNOLOGY", "YEAR"]).sum()


###
# private functions
###


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


def _format_ember_capacity_data(ember: pd.DataFrame) -> pd.DataFrame:
    """Formats data into otoole compatiable data structure"""

    df = ember.copy()

    df = df[(df.Category == "Capacity") & (df.Subcategory == "Fuel")].copy()
    df["TECHNOLOGY"] = df.Variable.map(TECHNOLOGY_MAPPER)
    df["TECHNOLOGY"] = df.TECHNOLOGY + df.COUNTRY
    df["REGION"] = "GLOBAL"
    df = df[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    return df.groupby(["REGION", "TECHNOLOGY", "YEAR"]).sum()


def _format_ember_generation_data(ember: pd.DataFrame) -> pd.DataFrame:
    """Formats data into otoole compatiable data structure"""

    df = ember.copy()

    df = df[
        (df.Category == "Electricity generation") & (df.Subcategory == "Fuel")
    ].copy()
    df["TECHNOLOGY"] = df.Variable.map(TECHNOLOGY_MAPPER)
    df["TECHNOLOGY"] = df.TECHNOLOGY + df.COUNTRY
    df["REGION"] = "GLOBAL"

    # TWh -> PJ
    # 1 TWh * (1PWh / 1000TWh) * (3600sec / hr) = 1 PWs = 1 PJ
    df["VALUE"] = df.VALUE.mul(3.6)
    df = df[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    return df.groupby(["REGION", "TECHNOLOGY", "YEAR"]).sum()
