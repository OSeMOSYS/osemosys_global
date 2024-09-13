"""Functions to perform demand regression"""

from typing import Optional
from sklearn.linear_model import LinearRegression
import pandas as pd
from spatial import get_spatial_mapping_country
from data import (
    get_historical_gdp_ppp_wb,
    get_historical_urban_pop_wb,
    get_historical_ember_demand,
)
from constants import SPATIAL_RESOLUTION


def _create_regression_dataframe(
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


def perform_regression(
    plexos: pd.DataFrame, ember: pd.DataFrame, urban: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """Performs the Linear Regression"""

    df = _create_regression_dataframe(plexos, ember, urban)

    urbanization = True if isinstance(urban, pd.DataFrame) else False

    # Groups the entries by <Spatial_Resolution> and calculates the regional linear
    # fit based on all historical values

    dfs = []

    for country in df.index.unique():

        df_country = df.loc[country].copy()

        if urbanization:
            country_lr = _regression_with_urbanization(df_country)
        else:
            country_lr = _regression_without_urbanization(df_country)

        dfs.append(country_lr)

    return pd.concat(dfs)


def _regression_with_urbanization(df: pd.DataFrame) -> pd.DataFrame:
    """Perform regression with urbanization

    If Urbanization is included linear regression occurs with multiple independent
    variables (GDPppp and % Urban population) for the dependent variable (Electricity demand).
    """

    lr = LinearRegression()

    lr.fit(
        df[["WB_GDPppp", "WB_Urb"]],
        df["ember_Elec"],
    )

    (
        df["intercept"],
        df["coef_GDPppp"],
        df["coef_Urb"],
    ) = [
        lr.intercept_,
        float(lr.coef_[0]),
        float(lr.coef_[1]),
    ]

    df["R2_GDPppp_Urb/Elec"] = lr.score(
        df[["WB_GDPppp", "WB_Urb"]],
        df["ember_Elec"],
    )

    return df


def _regression_without_urbanization(df: pd.DataFrame) -> pd.DataFrame:
    """Perform regression without urbanization

    If Urbanization is not included linear regression occurs with single independent
    variables (GDPppp) for the dependent variable (Electricity demand).
    """

    lr = LinearRegression()

    lr.fit(df[["WB_GDPppp"]], df["ember_Elec"])

    df["intercept"], df["coef_GDPppp"] = [
        lr.intercept_,
        float(lr.coef_),
    ]

    df["R2_GDPppp/Elec"] = lr.score(df[["WB_GDPppp"]], df["ember_Elec"])

    return df
