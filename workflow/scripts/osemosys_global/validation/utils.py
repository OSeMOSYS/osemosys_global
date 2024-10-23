"""Helper functions"""

import pandas as pd
from typing import Optional
import matplotlib.pyplot as plt
import datetime


def _join_data(
    modelled: pd.DataFrame, actual: pd.DataFrame, dataset_name: Optional[str] = None
) -> pd.DataFrame:
    """Joins modelled and actual data on indices"""

    if not dataset_name:
        dataset_name = "ACTUAL"

    modelled = modelled.rename(columns={"VALUE": "OSeMOSYS"})
    actual = actual.rename(columns={"VALUE": dataset_name})
    df = modelled.join(actual)

    assert len(df.index.get_level_values("REGION").unique()) == 1

    return df.droplevel("REGION")


def _get_last_validation_year() -> int:
    """Gets year to validate data up to.

    Note, this does not guarantee data will be avaialble from the dataset for this year.
    """
    # minus 1 as we can not validate current year
    return datetime.date.today().year - 1


###
# plotters
###


def plot_gen_cap(
    modelled: pd.DataFrame,
    actual: pd.DataFrame,
    variable: str,
    dataset_name: Optional[str] = None,
) -> dict[str, tuple[plt.figure, plt.axes]]:
    """Plots generation and capacity data"""

    assert modelled.index.names == actual.index.names

    if variable == "generation":
        units = "PJ"
    elif variable == "capacity":
        units = "GW"
    else:
        raise ValueError(
            f"Variable must be one of ['generation', 'capacity']. Recieved {variable}"
        )

    df = _join_data(modelled, actual, dataset_name).reset_index()
    df["TECH"] = df["TECHNOLOGY"].str[0:3]
    df["COUNTRY"] = df["TECHNOLOGY"].str[3:]

    data = {}

    countries = df.COUNTRY.unique()
    for country in countries:
        df_country = df[df.COUNTRY == country]
        years = df_country.YEAR.unique()
        n_rows = len(years)
        fig, axs = plt.subplots(n_rows, 1, figsize=(10, n_rows * 4))
        for i, year in enumerate(years):
            df_year = (
                df_country[df_country.YEAR == year]
                .drop(columns=["TECHNOLOGY", "YEAR", "COUNTRY"])
                .set_index("TECH")
            )
            title = f"{country} {variable.capitalize()} in {year}"
            if n_rows > 1:
                ax = axs[i]
            else:
                ax = axs
            df_year.plot(
                kind="bar", ax=ax, rot=45, title=title, xlabel="", ylabel=units
            )

        data[country] = (fig, axs)

    return data


def plot_emissions(
    modelled: pd.DataFrame,
    actual: pd.DataFrame,
    variable: str,
    dataset_name: Optional[str] = None,
) -> dict[str, tuple[plt.figure, plt.axes]]:
    """Plots generation and capacity data"""

    assert modelled.index.names == actual.index.names

    if variable == "emissions":
        units = "MT"
    elif variable == "emission_intensity":
        units = "g/kWh"
    else:
        raise ValueError(
            f"Variable must be one of ['emissions', 'emission_intensity']. Recieved {variable}"
        )

    df = _join_data(modelled, actual, dataset_name).reset_index()

    data = {}

    countries = df.EMISSION.unique()  # emission column holds country
    for country in countries:
        df_country = (
            df[df.EMISSION == country].drop(columns=["EMISSION"]).set_index("YEAR")
        )
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        title = f"{country} {variable.capitalize()}"
        df_country.plot(kind="bar", ax=ax, rot=45, title=title, xlabel="", ylabel=units)

        data[country] = (fig, ax)

    return data


###
# formatters
###


def format_rty_results(og: pd.DataFrame, mapper: dict[str, str]) -> pd.DataFrame:
    """Formats OG results for comparison on region, tech, year

    Mapper is to group different technologies together, to match external
    dataset aggregation.

    Works on:
    - ProductionByTechnologyAnnual
    - TotalCapacityAnnual
    """

    df = og.copy()

    year = _get_last_validation_year()

    if len(df.columns) == 1:
        df = df.reset_index()

    df = df[(df.TECHNOLOGY.str.startswith("PWR")) & (df.YEAR <= year)]
    df["COUNTRY"] = df.TECHNOLOGY.str[6:9]
    df["CODE"] = df.TECHNOLOGY.str[3:6]
    df["CODE"] = df.CODE.map(mapper)
    df = df.dropna(subset="CODE")
    df["TECHNOLOGY"] = df.CODE + df.COUNTRY
    df = df[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    return df.groupby(["REGION", "TECHNOLOGY", "YEAR"]).sum()


def format_rey_results(og: pd.DataFrame) -> pd.DataFrame:
    """Formats OG results for comparison on region, emission, year

    Works on:
    - AnnualEmissions
    """

    df = og.copy()

    year = _get_last_validation_year()

    if len(df.columns) == 1:
        df = df.reset_index()

    # emission is used to track country
    df = df[df.YEAR <= year]
    df["EMISSION"] = df.EMISSION.str[3:6]
    return df.set_index(["REGION", "EMISSION", "YEAR"])
