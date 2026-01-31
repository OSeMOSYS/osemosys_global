"""Module for reading in data sources"""

import pandas as pd


def import_plexos_2015(f: str) -> pd.DataFrame:
    """Imports PLEXOS-World 2015 model file as basis for the spatial mapping.
    
    PLEXOS_World_2015_Gold_V1.1.xlsx
    """
    return pd.read_excel(f, sheet_name="Memberships")


def import_iamc(f: str) -> pd.DataFrame:
    """Imports SSP GDPppp and Population projections from IIASA
    https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=30
    
    iamc_db_GDPppp_Countries.xlsx
    iamc_db_POP_Countries.xlsx
    iamc_db_URB_Countries.xlsx
    """
    return pd.read_excel(f)


def import_iamc_missing(f: str, metric: str) -> dict[str, pd.DataFrame]:
    """Imports custom GDPppp and Population projections for countries not included in SSP datasets
    
    iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx
    """

    if metric.lower() == "gdp":
        sheet_name = "GDP|PPP"
    elif metric.lower() == "pop":
        sheet_name = "POP"
    elif metric.lower() == "urb":
        sheet_name = "URB"
    else:
        raise NotImplementedError

    return pd.read_excel(f, sheet_name=sheet_name).set_index("Region")


def import_td_losses(f: str) -> pd.DataFrame:
    """Imports transmission and distribution losses projections

    https://www.sciencedirect.com/science/article/pii/S0142061518335075?via%3Dihub
    
    T&D Losses.xlsx
    """
    return pd.read_excel(f)


def import_hourly_demand(f: str) -> pd.DataFrame:
    """Imports PLEXOS-World 2015 hourly demand data
    
    All_Demand_UTC_2015.csv
    """

    return pd.read_csv(f, encoding="latin-1")


def import_ember_elec(f: str) -> pd.DataFrame:
    """Imports EMBER yearly electricity data
    
    ember_yearly_electricity_data.csv
    """

    df = pd.read_csv(f, encoding="latin-1")
    return (
        df[["ISO 3 code", "Year", "Variable", "Value"]]
        .rename(columns={"Value": "ember_Elec", "ISO 3 code": "Country"})
        .dropna()
        .set_index("Country")
    )

