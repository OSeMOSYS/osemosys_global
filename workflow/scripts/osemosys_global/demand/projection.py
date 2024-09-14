"""Functions for perfoming the demand projections"""

import pandas as pd
from data import get_iamc_data
from regression import get_regression_coefficients
from constants import START_YEAR, END_YEAR

from typing import Optional


def _get_year_interval(start: int, end: int) -> range:
    """Standard interval for demand projections"""
    return range(start, end + 1, 5)


def perform_country_projection(
    plexos: pd.DataFrame,
    iamc_gdp: pd.DataFrame,
    iamc_pop: pd.DataFrame,
    iamc_urb: pd.DataFrame,
    iamc_missing: pd.DataFrame,
    lr: pd.DataFrame,
) -> pd.DataFrame:
    """Perfroms demand projections at country level

    Projects electricity demand by making use of historic relationships and SSP
    specific Population, GDP|PPP and optionally urbanization projections
    """

    gdp = get_iamc_data(plexos, iamc_gdp, iamc_missing, "gdp")
    pop = get_iamc_data(plexos, iamc_pop, iamc_missing, "pop")
    urb = get_iamc_data(plexos, iamc_urb, iamc_missing, "urb")

    lr_coef = get_regression_coefficients(lr, True)

    base = _get_base_data(gdp, pop, lr_coef)
    
    df = get_electrical_projection_country(base, urb)
    
    return df

def _get_base_data(
    iamc_gdp: pd.DataFrame,
    iamc_pop: pd.DataFrame,
    lr_coef: pd.DataFrame,
) -> pd.DataFrame:
    """Creates base dataframe for demand projections with the required coefficients"""

    base = iamc_gdp[["child_object", "Scenario"]]

    df = pd.merge(
        base,
        lr_coef,
        left_on="child_object",
        right_index=True,
        how="left",
    )

    # Divides the country level GDP|PPP data (converted from billions to millions)
    # with the population (millions) to get GDP|PPP pp.

    years = _get_year_interval(START_YEAR, END_YEAR)

    for year in years:
        df[year] = (base[year] * 1000) / iamc_pop[year]

    return df


def get_electrical_projection_country(
    base: pd.DataFrame, iamc_urb: Optional[pd.DataFrame] = None
) -> pd.DataFrame:
    """Country-level GDP|PPP pp values to project country-level electricity demand pp"""

    if isinstance(iamc_urb, pd.DataFrame):
        assert "coef_Urb" in base.columns
        cols = ["child_object", "Scenario", "coef_GDPppp", "coef_Urb", "intercept"]
    else:
        cols = ["child_object", "Scenario", "coef_GDPppp", "intercept"]

    df = base[cols].copy()

    df["Variable"] = "Demand|projected|pp"

    years = _get_year_interval(START_YEAR, END_YEAR)

    if isinstance(iamc_urb, pd.DataFrame):
        for year in years:
            df[year] = (
                df["coef_GDPppp"] * base[year]
                + df["coef_Urb"] * iamc_urb[year]
                + df["intercept"]
            )
    else:
        for year in years:
            df[year] = df["coef_GDPppp"] * base[year] + df["intercept"]

    return df
