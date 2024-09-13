"""Functions to extract relevent data"""

import pandas as pd
import wbgapi as wb
from datetime import datetime

from spatial import get_spatial_mapping_country
from constants import (
    POP_COUNTRIES_SOURCE,
    GDP_PPP_COUNTRIES_SOURCE,
    URB_COUNTRIES_SOURCE,
    PATHWAY,
    SPATIAL_RESOLUTION,
)


def get_nodal_plexos_demand(plexos: pd.DataFrame) -> pd.DataFrame:
    """Gets historical nodal demand

    Determines relative 2015 share of demand per sub-country node
    """

    raw = plexos.copy()

    # Sums the hourly demand as retrieved from the PLEXOS-World dataset to year total (in MWh) and drops all hourly values.
    df = raw.drop(columns=["Datetime"])
    df.loc["Node_Demand_2015"] = df.sum()
    df = df.iloc[8760:]

    # Transposes the dataframe and uses the original headers as column entry.
    df = df.transpose().reset_index().rename(columns={"index": "PLEXOS_Nodes"})

    # Adds country entry to dataframe (e.g. NA-USA-CA Node gets column entry NA-USA country)
    df.insert(
        loc=0,
        column="Country",
        value=df.PLEXOS_Nodes.str.split("-", expand=True)[1],
    )

    df.insert(loc=1, column="PLEXOS_Countries", value=df["PLEXOS_Nodes"].str[:6])

    # Creates a dataframe excluding all sub-country nodes
    country = (
        df[df["PLEXOS_Countries"] == df["PLEXOS_Nodes"]]
        .drop(columns=["PLEXOS_Nodes"])
        .rename(columns={"Node_Demand_2015": "Country_Demand_2015"})
    )

    # Adds country-level 2015 demand in column adjacent to sub-country level 2015 demand
    # and calculates relative share per country.

    nodes = pd.merge(
        df,
        country[["PLEXOS_Countries", "Country_Demand_2015"]],
        on="PLEXOS_Countries",
        how="left",
    )

    nodes["Share_%_Country_Demand"] = (
        nodes["Node_Demand_2015"] / nodes["Country_Demand_2015"]
    )

    return nodes


def get_historical_gdp_ppp_wb(long: bool = True) -> pd.DataFrame:
    """Gets historical GDPppp per capita from the World Bank API"""
    df = _extract_wb("NY.GDP.PCAP.PP.KD")
    if long:
        return _longify_wb(df, "WB_GDPppp")
    else:
        return df


def get_historical_urban_pop_wb(long: bool = True) -> pd.DataFrame:
    """Gets historical Urban population (% of total population) from the World Bank API"""
    df = _extract_wb("SP.URB.TOTL.IN.ZS")
    if long:
        return _longify_wb(df, "WB_Urb")
    else:
        return df


def _extract_wb(column: str) -> pd.DataFrame:
    """Extracts world bank data

    2000 used as first year as EMBER dataset doesn't go back further.
    """
    return (
        wb.data.DataFrame([column], mrv=datetime.now().year - 2000)
        .reset_index()
        .rename(columns={"economy": "Country"})
        .set_index("Country")
    )


def _longify_wb(df: pd.DataFrame, column_name: str) -> pd.DataFrame:
    """Converts world bank data into long format

    TODO: Change this to a melt function
    """

    dfs = []

    for year in df.columns:
        data = df[[year]].rename(columns={year: column_name})
        data["Year"] = year.replace("YR", "")
        dfs.append(data)

    return pd.concat(dfs)


def get_historical_ember_demand(ember: pd.DataFrame) -> pd.DataFrame:
    """Gets historical ember electricity data per capita"""

    df = ember.copy()

    # Electricity demand per capita only
    df = df[df.Variable == "Demand per capita"].drop(columns={"Variable"})

    # Conversion to kWh
    df["ember_Elec"] = df["ember_Elec"] * 1000
    df["Year"] = df["Year"].astype("str")  #################### why

    return df


def _iamc_data_available(
    iamc: pd.DataFrame, iamc_missing: pd.DataFrame, spatial_mapping: pd.DataFrame
) -> bool:
    """Checks whether pathway data is available for all included countries.

    In case pathway data is not available for all countriesit checks wheter custom data is provided in 'iamc_db_POP_POP_Countries_Missing.xlsx'. If not, an error POPs up indicatingfor which countries data is missing. Data has to be manually added to project demand for all countries.
    """

    country_missing = spatial_mapping[(~spatial_mapping.index.isin(iamc.index))]

    for x in country_missing.index:
        if not x in iamc_missing.index:
            return False
    return True


def get_iamc_data(
    plexos: pd.DataFrame, iamc: pd.DataFrame, iamc_missing: pd.DataFrame, metric: str
) -> pd.DataFrame:
    """Gets full iamc data"""

    if metric == "gdp":
        model = GDP_PPP_COUNTRIES_SOURCE
    elif metric == "pop":
        model = POP_COUNTRIES_SOURCE
    elif metric == "urb":
        model = URB_COUNTRIES_SOURCE
    else:
        raise NotImplementedError

    spatial_mapping = get_spatial_mapping_country(plexos)

    df_original = iamc[
        (iamc["Model"] == model) & (iamc["Scenario"] == PATHWAY)
    ].set_index("Region")

    if not _iamc_data_available(df, iamc_missing, spatial_mapping):
        raise ValueError("Country data for is not available in custom dataset!")

    df_missing = iamc_missing[(iamc_missing["Scenario"] == PATHWAY)]

    df = pd.concat([df_original, df_missing])

    # Filters data for relevant to be modelled countries
    return pd.merge(
        spatial_mapping[[SPATIAL_RESOLUTION]],
        df,
        left_index=True,
        right_index=True,
        how="inner",
    )
