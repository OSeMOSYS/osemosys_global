"""Functions for perfoming the demand projections"""

import pandas as pd
import numpy as np
from regression import get_regression_coefficients
from constants import START_YEAR, END_YEAR, PEAK_RATIO_FACTOR
from spatial import get_spatial_mapping_country, get_spatial_mapping_node
from data import get_nodal_plexos_demand

from typing import Optional

###
# Public functions
###


def perform_node_projections(
    lr: pd.DataFrame,
    plexos: pd.DataFrame,
    plexos_demand: pd.DataFrame,
    iamc_gdp: pd.DataFrame,
    iamc_pop: pd.DataFrame,
    iamc_urb: pd.DataFrame,
    td_losses: pd.DataFrame,
) -> pd.DataFrame:
    """Creates yearly demand projections for all nodes"""

    df = perform_node_projection_step(
        lr, plexos, plexos_demand, iamc_gdp, iamc_pop, iamc_urb, td_losses
    )

    proj = _interpolate_yearly_demand(df)
    proj = _get_node_peak_demand_ratio(plexos_demand, proj)
    # return _adjust_with_peak_demand(proj, PEAK_RATIO_FACTOR)
    return proj


def perform_node_projection_step(
    lr: pd.DataFrame,
    plexos: pd.DataFrame,
    plexos_demand: pd.DataFrame,
    iamc_gdp: pd.DataFrame,
    iamc_pop: pd.DataFrame,
    iamc_urb: pd.DataFrame,
    td_losses: pd.DataFrame,
) -> pd.DataFrame:
    """Performs demand projections at NODE level in step intervals (5 years)

    TODO: Would be nice to include the td_losses as an optional argument
    """

    df = perform_country_projection_step(lr, iamc_gdp, iamc_pop, iamc_urb)

    plexos_node_demand = get_nodal_plexos_demand(plexos_demand)

    df = _apply_td_losses(plexos, plexos_node_demand, td_losses, df)

    return _downscale_demand(plexos, plexos_node_demand, df)


def perform_country_projection_step(
    lr: pd.DataFrame,
    iamc_gdp: pd.DataFrame,
    iamc_pop: pd.DataFrame,
    iamc_urb: pd.DataFrame,
    per_capita: Optional[bool] = False,
    **kwargs
) -> pd.DataFrame:
    """Performs demand projections at COUNTRY level in step intervals (5 years)

    Projects electricity demand by making use of historic relationships and SSP
    specific Population, GDP|PPP and optionally urbanization projections
    """

    lr_coef = get_regression_coefficients(lr, True)

    base = _get_base_data(iamc_gdp, iamc_pop, lr_coef)

    if per_capita:
        df = _get_electrical_projection_country(base, iamc_urb)
    else:
        df = _get_electrical_projection_country(base, iamc_urb, iamc_pop)

    return df


###
# Private functions
###


def _get_year_interval(start: int, end: int) -> range:
    """Standard interval for demand projections"""
    return range(start, end + 5, 5)


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
        df[year] = (iamc_gdp[year] * 1000) / iamc_pop[year]

    return df


def _get_electrical_projection_country(
    base: pd.DataFrame,
    iamc_urb: Optional[pd.DataFrame] = None,
    iamc_pop: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Country-level GDP|PPP pp values to project country-level electricity demand pp

    Data is returned as a per-capita value
    """

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

    if isinstance(iamc_pop, pd.DataFrame):
        return _convert_per_capita_to_total(df, iamc_pop)
    else:
        return df


def _convert_per_capita_to_total(
    df_per_capita: pd.DataFrame, iamc_pop: pd.DataFrame
) -> pd.DataFrame:
    """Converts the per-capita projection to total demand

    Multiplies the country-level projected demand pp (in kWh) with the total
    population (in millions) to get country-level total projected demand (in GWh).
    """

    df = pd.DataFrame()
    years = _get_year_interval(START_YEAR, END_YEAR)
    for year in years:
        df[year] = df_per_capita[year] * iamc_pop[year]

    return df


def _apply_td_losses(
    plexos: pd.DataFrame,
    plexos_node_demand: pd.DataFrame,
    td_losses: pd.DataFrame,
    projection: pd.DataFrame,
) -> pd.DataFrame:
    """Apply Transmission Distribution losses to country level demand

    Explicit modelling of domestic transmission and distribution is not incorporated
    in PLEXOS-World. Country-level T&D losses per 5-year interval are added to the
    projected electricity demand based on
    Sadovskaia et al., 2019; https://doi.org/10.1016/j.ijepes.2018.11.012.
    Study includes data for up till 2050. 5-year intervals after that are manually
    added (Maarten Brinkerink) and values kept equal compared to 2050.
    """

    df = td_losses.copy()

    spatial_mapping = get_spatial_mapping_country(plexos)

    # Checks whether T&D is available for all included countries
    losses_missing = spatial_mapping[(~spatial_mapping.index.isin(df.Country))]

    assert losses_missing.empty

    df = df.set_index("Country")

    ctry_losses = df[df.index.isin(spatial_mapping.index)].reindex(
        spatial_mapping.index
    )

    years = _get_year_interval(START_YEAR, END_YEAR)

    projection_with_losses = pd.DataFrame()

    # Add T&D losses to the projected country-level demand.
    for year in years:
        projection_with_losses[year] = (
            (projection[year] * ctry_losses[year] / 100 + projection[year])
            .round(2)
            .astype(float)
        )

    # Constraints the forecasted final demand to 2015 baseline values as minimum
    # In case of linear regression, smaller countries with signficantly lower
    # projected independent variables (GDP, Urbanization) compared to the regional
    # average can lead to very low and often negative projected demand
    # values (e.g. EU-KOS). Hence, a comparison is being made to the 2015 baseline
    # demand values with the assumption that a decline in electricity demand is
    # not realistic (note: as of now no decoupling of GDP growth and energy
    # demand reduction has been assumed).

    dfs = []

    for country in projection_with_losses.index.unique():
        ctry_demand = projection_with_losses.loc[country]
        node_demand = plexos_node_demand[
            plexos_node_demand.Country.str.contains(country)
        ]

        total = ctry_demand.clip((node_demand.iloc[0]["Country_Demand_2015"] / 1000))

        dfs.append(pd.DataFrame(total).T)

    return pd.concat(dfs)


def _downscale_demand(
    plexos: pd.DataFrame, plexos_node_demand: pd.DataFrame, projection: pd.DataFrame
) -> pd.DataFrame:
    """Downscales country demand to nodal demand"""

    df = projection.copy()

    nodes = plexos_node_demand[
        ["Country", "PLEXOS_Countries", "PLEXOS_Nodes", "Share_%_Country_Demand"]
    ].copy()

    df = pd.merge(
        nodes,
        df,
        left_on="Country",
        right_index=True,
    )

    # Downscales projected country level demand by the 2015 shares of sub-country nodes as proxy.

    years = _get_year_interval(START_YEAR, END_YEAR)
    for year in years:
        df[year] = (df[year] * df["Share_%_Country_Demand"]).round(2).astype(float)

    df = df.set_index("PLEXOS_Nodes")

    spatial_mapping = get_spatial_mapping_node(plexos)

    # Filters data for the to-be modelled counties.
    df = df[(df.index.isin(spatial_mapping.index))]

    df.insert(loc=3, column="Unit", value="GWh")

    return df


def _interpolate_yearly_demand(projection: pd.DataFrame) -> pd.DataFrame:
    """Interpolates data to yearly datapoints

    Interpolates 5-yearly values to yearly values and determines final electricity demand
    per node per year
    """

    df = projection.copy()

    step_years = _get_year_interval(START_YEAR, END_YEAR)
    all_years = range(START_YEAR, END_YEAR + 1)

    for year in all_years:
        if year not in step_years:
            df[year] = np.NaN

    # Filters columns with numerical dtypes.
    df_interp = df[all_years].copy()

    # Linearly interpolates values between years.
    df_interp = df_interp.interpolate(method="linear", axis=1)

    # Merges dataframes.

    return pd.merge(
        df[["Country", "PLEXOS_Countries", "Share_%_Country_Demand", "Unit"]],
        df_interp,
        left_index=True,
        right_index=True,
    )


def _get_node_peak_demand_ratio(
    plexos_demand: pd.DataFrame, projection: pd.DataFrame
) -> pd.DataFrame:

    demand = plexos_demand.copy().drop(columns=["Datetime"])
    demand = pd.DataFrame(demand.max()).reset_index()
    demand.columns = ["PLEXOS_Nodes", "Node_Peak_Demand_2015"]

    plexos_demand_node = get_nodal_plexos_demand(plexos_demand)

    # Calculates 2015 peak demand / total demand ratio
    demand["Ratio_Peak/Total_Demand_2015"] = (
        demand["Node_Peak_Demand_2015"] / plexos_demand_node["Node_Demand_2015"]
    )

    demand = demand[["PLEXOS_Nodes", "Ratio_Peak/Total_Demand_2015"]]

    # Merges dataframes.
    return (
        pd.merge(
            demand,
            projection,
            left_on="PLEXOS_Nodes",
            right_index=True,
            how="right",
        )
        .set_index("PLEXOS_Nodes")
        .drop(columns={"Share_%_Country_Demand"})
    )


def _adjust_with_peak_demand(
    projection: pd.DataFrame, peak_ratio_factor: float
) -> pd.DataFrame:
    """Calculates projected hourly peak demand by using the relative 2015 peak demand as proxy.

    The peak to total demand ratio can be adjusted by changing the peak_ratio_factor.
    """

    df = projection.copy()

    for year in range(START_YEAR, 2101):  # TODO fix this to END_YEAR
        df[year] = (
            df[year] * peak_ratio_factor * df["Ratio_Peak/Total_Demand_2015"] * 1000
        ).round(2)

    df["Unit"] = "MW"

    return df
