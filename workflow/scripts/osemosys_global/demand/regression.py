"""Functions to perform demand regression"""

from typing import Optional

import pandas as pd
from .spatial import get_spatial_mapping_country
from .data import (
    get_historical_gdp_ppp_wb,
    get_historical_urban_pop_wb,
    get_historical_ember_demand,
)
from .constants import SPATIAL_RESOLUTION


def create_regression_dataframe(
    plexos: pd.DataFrame, ember: pd.DataFrame, urban: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """
    Creates dataframe to be used in the linear regression
    """

    gdp_wb = get_historical_gdp_ppp_wb(long=True)[["Year", "WB_GDPppp"]]
    dem_ember = get_historical_ember_demand(ember)[["Year", "ember_Elec"]]

    df = pd.merge(
        gdp_wb,
        dem_ember,
        left_on=["Year", "Country"],
        right_on=["Year", "Country"],
    )

    if isinstance(urban, pd.DataFrame):
        urb_wb = get_historical_urban_pop_wb(long=True)[["Year", "WB_Urb"]]
        df = pd.merge(
            df,
            urb_wb,
            left_on=["Year", "Country"],
            right_on=["Year", "Country"],
        )

    # Drops all entries that don't have an inner match
    df = df.dropna()

    spatial_mapping = get_spatial_mapping_country(plexos)[
        ["parent_object", "child_object"]
    ]

    df = pd.merge(
        df,
        spatial_mapping,
        left_index=True,
        right_index=True,
    )

    return df.set_index(SPATIAL_RESOLUTION)


