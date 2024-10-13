"""Processes data from IRENA IRENASTAT database

https://www.irena.org/Data/Downloads/IRENASTAT
"""

import pandas as pd
from typing import Optional

###
# constants for data allignment
###

TECHNOLOGY_MAPPER = {
    "Solar photovoltaic": "SPV",
    "Onshore wind energy": "WON",
    "Offshore wind energy": "WOF",
    "Renewable hydropower": "HYD",
    "Mixed Hydro Plants": "HYD",
    "Marine energy": "WAV",
    "Solid biofuels": "BIO",
    "Renewable municipal waste": "WAS",
    "Liquid biofuels": "BIO",
    "Biogas": "BIO",
    "Geothermal energy": "GEO",
    "Pumped storage": "PHP",
    "Coal and peat": "COA",
    "Oil": "OIL",
    "Natural gas": "GAS",
    "Nuclear": "URN",
}

###
# public functions
###


def get_irena_capacity(csv_file: str, iso_codes: str, **kwargs) -> pd.DataFrame:
    df = _read_irena_data(csv_file, iso_codes)
    return _format_irena_capacity_data(df)


def get_irena_generation(csv_file: str, iso_codes: str, **kwargs) -> pd.DataFrame:
    df = _read_irena_data(csv_file, iso_codes)
    return _format_irena_generation_data(df)


def format_og_generation(prod_tech_annual: pd.DataFrame) -> pd.DataFrame:
    """Formats ProductionByTechnologyAnnual data for irena comparison"""

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
        "WON": "WON",
        "WOF": "WOF",
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


def _read_irena_data(csv_file: str, iso_codes: Optional[str] = None) -> pd.DataFrame:
    """Reads *.csv IRENA data from https://www.irena.org/Data/Downloads/IRENASTAT

    Power Capacity and Generation by Country/area, TEchnology
    """

    def _get_iso_mapper(iso_codes: str) -> dict[str, str]:
        df = pd.read_csv(iso_codes)
        return df.set_index("name")["alpha-3"].to_dict()

    converter = lambda x: x.replace("-", "0")
    df = pd.read_csv(
        csv_file,
        skiprows=2,
        converters={"Electricity statistics (MW/GWh)": converter},
        encoding="latin-1",
    )
    df = df.rename(
        columns={
            "Country/area": "COUNTRY",
            "Technology": "TECHNOLOGY",
            "Year": "YEAR",
            "Electricity statistics (MW/GWh)": "VALUE",
        }
    )
    df["VALUE"] = df.VALUE.astype(float)

    if iso_codes:
        iso_map = _get_iso_mapper(iso_codes)
        df["COUNTRY"] = df.COUNTRY.map(iso_map)

    df["TECHNOLOGY"] = df.TECHNOLOGY.map(TECHNOLOGY_MAPPER)

    return df


def _format_irena_capacity_data(eia: pd.DataFrame) -> pd.DataFrame:
    """Formats data into otoole compatiable data structure"""
    raise NotImplementedError


def _format_irena_generation_data(irena: pd.DataFrame) -> pd.DataFrame:
    """Formats data into otoole compatiable data structure"""

    df = irena.copy()

    assert len(df["Data Type"].unique()) == 1

    df["TECHNOLOGY"] = df.TECHNOLOGY + df.COUNTRY
    df["REGION"] = "GLOBAL"

    # GWh -> PJ
    # 1 GWh * (1TWh / 1000GWh) * (1PWh / 1000TWh) * (3600sec / hr) = 1 PWs = 1 PJ
    df["VALUE"] = df.VALUE.div(1000000).mul(3600)
    df = df[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    return df.groupby(["REGION", "TECHNOLOGY", "YEAR"]).sum()
