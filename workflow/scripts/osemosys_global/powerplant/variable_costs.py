"""Generates variable costs"""

import logging

logger = logging.getLogger(__name__)

import pandas as pd
from typing import Optional

from read import (
    import_cmo_forecasts,
    import_fuel_prices,
    import_set,
    import_fuel_limits,
)

# can not use constant file defintions due to waste being tagged incorrectly
MINING_FUELS = ["COA", "COG", "GAS", "PET", "URN", "OIL", "OTH", "WAS"]
RENEWABLE_FUELS = ["BIO", "GEO", "HYD", "SPV", "CSP", "WAV", "WON", "WOF"]

from constants import (
    CMO_DATA_YEAR,
    BIOMASS_VAR_COSTS,
    NUCLEAR_VAR_COSTS,
    WASTE_VAR_COSTS,
    INT_COST_FACTOR,
)


def calculate_cmo_forecasts(cmo: pd.DataFrame, data_year: int) -> pd.DataFrame:
    """
    Data for variable costs of fuels taken from:

    World Bank Commodity Market Outlooks
    https://www.worldbank.org/en/research/commodity-markets

    ORIGINAL UNITS:
           MINCOA MINOIL MINGAS  KFNGAS_US KFNGAS_JP
    Unit   $/mt   $/bbl  $/mmbtu $/mmbtu   $/mmbtu
    """

    df = (
        cmo[["Commodity", data_year]]
        .replace(
            {
                "Coal, Australia": "COA",
                "Crude oil, Brent": "OIL",
                "Natural gas, Europe": "KFNGAS_EU",
                "Natural gas, U.S.": "KFNGAS_US",
                "Liquefied natural gas, Japan": "KFNGAS_JP",
            },
            regex=True,
        )
        .rename(columns={"Commodity": "FUEL", data_year: "VALUE"})
    )

    # Add average GAS price based on the 'KFNGAS_EU', KFNGAS_US' and KFNGAS_JP' entries
    gas_avg = df.loc[df["FUEL"].str.contains("GAS")].loc[:, "VALUE"].mean()
    df.loc[len(df.index)] = ["GAS", gas_avg]
    df = df.loc[~df["FUEL"].isin(["KFNGAS_EU", "KFNGAS_US", "KFNGAS_JP"])]

    """
    Convert to $/PJ (Values taken from kylesconverter.com):
      1 mt coal contains 29.31 GJ. We want $mill/PJ so divide by 29.31
      1 mmbtu = 0.000001055056 PJ. We want $mill/PJ so divide by 0.000001055056 
          and multiply with 1000000 (= / 1.055056).
      1 bbl = 0.00000612 PJ (Barrels of Oil) Original was $/bbl so divide by 
          0.00000612 and divide by 1000000 ( = / 6.12)
    """
    energy_mapper = {
        "GAS": 1.055056,
        "COA": 29.31,
        "OIL": 6.12,
    }

    unit_mapper = {
        "GAS": "$/mmbtu",
        "COA": "$/mt",
        "OIL": "$/bbl",
    }

    df["COUNTRY"] = "INT"
    df["UNIT"] = df.FUEL.map(unit_mapper)
    df["ENERGY_CONTENT"] = df.FUEL.map(energy_mapper)

    return df


def copy_cmo_fuel_price(
    cmo: pd.DataFrame, old_fuel: str, new_fuel: str
) -> pd.DataFrame:
    """Copies a CMO fuel and applies the price to a new fuel"""
    df = cmo.copy()
    return df.loc[df.FUEL == old_fuel].replace(old_fuel, new_fuel)


def get_static_fuel_price(
    fuel: str, country: str, unit: str, energy_content: int | float, value: int | float
) -> pd.DataFrame:
    """Gets static user defined fuel price"""

    return pd.DataFrame(
        [[fuel, country, unit, energy_content, value]],
        columns=["FUEL", "COUNTRY", "UNIT", "ENERGY_CONTENT", "VALUE"],
    )


def expand_cmo_data(cmo: pd.DataFrame, years: list[int]) -> pd.DataFrame:
    """Expands CMO data to match user defined years"""
    df = cmo.copy()

    for year in years:
        df[str(year)] = df["VALUE"]

    return df.drop(columns="VALUE").reset_index(drop=True)


def get_user_fuel_years(user_fuel_prices: pd.DataFrame) -> list[int]:
    """Gets years that user fuel prices are defined over"""

    df = user_fuel_prices.copy()
    df = df.drop(columns=["FUEL", "COUNTRY", "UNIT", "ENERGY_CONTENT"])

    try:
        return [int(x) for x in df.columns]
    except ValueError as ex:
        logger.error("Unexpected column header in user defined fuel prices.")
        raise ValueError(ex)


def merge_cmo_user(cmo: pd.DataFrame, user: pd.DataFrame) -> pd.DataFrame:
    """Merges CMO data with user defined data.

    User defined data is used wherever possible, with CMO being used to fill
    missing data
    """

    assert len(cmo.columns) == len(user.columns)

    user_copy = user.copy()
    cmo_copy = cmo[user_copy.columns]

    assert user_copy.columns.equals(cmo_copy.columns)

    df = pd.concat([user_copy, cmo_copy])

    df = df.drop_duplicates(subset=["FUEL", "COUNTRY"], keep="first")

    return df


def apply_energy_content(costs: pd.DataFrame) -> pd.DataFrame:
    """Applies energy content calculation to convert everything into M$/PJ"""

    df = costs.copy()

    # cols will be years
    cols = [
        x for x in df.columns if x not in ["FUEL", "COUNTRY", "ENERGY_CONTENT", "UNIT"]
    ]

    df = df.drop(columns=["UNIT"])

    df = pd.melt(
        df,
        id_vars=["FUEL", "COUNTRY", "ENERGY_CONTENT"],
        value_vars=cols,
        var_name="YEAR",
        value_name="VALUE",
    )

    df["VALUE"] = df.VALUE.div(df.ENERGY_CONTENT)

    df = df.drop(columns=["ENERGY_CONTENT"])

    return df


def apply_technology_name(costs: pd.DataFrame) -> pd.DataFrame:

    df = costs.copy()

    df.loc[df["FUEL"].isin(["BIO", "WAS"]), "TECHNOLOGY"] = (
        "RNW" + df["FUEL"] + df["COUNTRY"]
    )
    df.loc[~(df["FUEL"].isin(["BIO", "WAS"])), "TECHNOLOGY"] = (
        "MIN" + df["FUEL"] + df["COUNTRY"]
    )
    df["VALUE"] = df["VALUE"].round(2)

    return df


def get_countries_from_techs(technologies: pd.Series) -> list[str]:
    """Gets unique countries in the model"""

    techs = technologies.copy()
    techs = techs[techs.str.startswith("PWR")]
    techs = techs.map(lambda x: x[6:9])

    return techs.unique().tolist()


def get_nodes_from_techs(technologies: pd.Series) -> list[str]:
    """Gets unique nodes in the model"""

    techs = technologies.copy()
    techs = techs[techs.str.startswith("PWR")]
    techs = techs.map(lambda x: x[6:11])

    return techs.unique().tolist()


def filter_var_cost_technologies(
    df: pd.DataFrame, technologies: pd.Series
) -> pd.DataFrame:
    """Filters variable costs for available technologies

    In some cases, a powerplant may not exist in a country, and also not be investable
    (ie. Waste plants). In these situations, we want to filter them out
    """
    available_techs = technologies[technologies.str.startswith("MIN")]
    return df[df.index.get_level_values("TECHNOLOGY").isin(available_techs)].copy()


def assign_international_limits(
    var_costs: pd.DataFrame, fuel_limits: pd.DataFrame
) -> pd.DataFrame:
    """Checks if a mining tech can contribute to international markets

    If a country has a fuel limit, the country can contribute to international markets. If the
    country does not have a fuel limit, the international variable cost is set to a very high value.
    """

    df = var_costs.copy()

    cost = 999999

    # no international trading
    if fuel_limits.empty:
        df["VALUE"] = cost
        return df

    limits = fuel_limits.copy()
    limits["name"] = "MIN" + limits.FUEL + limits.COUNTRY
    limited_techs = limits.name.unique().tolist()

    df.loc[df.TECHNOLOGY.isin(limited_techs), "VALUE"] = cost

    return df


def _get_cmo_template(
    costs: pd.DataFrame,
    fuels: list[str],
    countries: list[str],
    years: list[int],
    suffix: Optional[str] = None,
) -> pd.DataFrame:
    """Gets variable costs defined by cmo template

    ie. This will parse the costs data for only international prices
    """

    # get year bounds as interpolation can use years outside model scope
    min_year = min(min(years), min(costs.YEAR.astype(int)))
    max_year = max(max(years), max(costs.YEAR.astype(int)))

    df = costs[costs.FUEL.isin(fuels)].copy()

    df = (
        df[df.COUNTRY == "INT"]
        .copy()
        .drop(columns="COUNTRY")
        .set_index(["FUEL", "YEAR"])
    )
    fuel = df.index.get_level_values("FUEL").unique().tolist()
    idx = pd.MultiIndex.from_product(
        [fuel, range(min_year, max_year + 1)], names=["FUEL", "YEAR"]
    )
    df = df.reindex(idx)
    df = df.interpolate(method="linear", limit_direction="both").reset_index()
    df["COUNTRY"] = [countries] * len(df)
    df = df.explode("COUNTRY")
    if suffix:
        df["TECHNOLOGY"] = suffix + df.FUEL + df.COUNTRY
    else:
        df["TECHNOLOGY"] = df.FUEL + df.COUNTRY

    # filter on years within actual model
    df = df[df.YEAR.isin(years)].copy()

    return df[["TECHNOLOGY", "YEAR", "VALUE"]]


def _get_user_template(
    costs: pd.DataFrame,
    fuels: list[str],
    years: list[int],
    suffix: Optional[str] = None,
) -> pd.DataFrame:
    """Gets variable costs defined by user

    ie. This will parse the costs data for only domestic prices
    """

    # get year bounds as interpolation can use years outside model scope
    min_year = min(min(years), min(costs.YEAR.astype(int)))
    max_year = max(max(years), max(costs.YEAR.astype(int)))

    df = costs[(costs.FUEL.isin(fuels)) & (costs.COUNTRY != "INT")].copy()

    if suffix:
        df["TECHNOLOGY"] = suffix + df.FUEL + df.COUNTRY
    else:
        df["TECHNOLOGY"] = df.FUEL + df.COUNTRY

    df = df.drop(columns=["COUNTRY", "FUEL"]).set_index(["TECHNOLOGY", "YEAR"])
    techs = df.index.get_level_values("TECHNOLOGY").unique().tolist()
    idx = pd.MultiIndex.from_product(
        [techs, range(min_year, max_year + 1)], names=["TECHNOLOGY", "YEAR"]
    )
    df = df.reindex(idx)
    df = df.interpolate(method="linear", limit_direction="both").reset_index()

    # filter on years within actual model
    df = df[df.YEAR.isin(years)].copy()

    return df[["TECHNOLOGY", "YEAR", "VALUE"]]


def get_domestic_mining_data(
    costs: pd.DataFrame,
    mining_fuels: list[str],
    countries: list[str],
    years: list[int],
    region: str,
) -> pd.DataFrame:
    """Gets domestic mining data (ie. operates on mode 1)

    Logic is:
    - Assigns data for all regions based on international data
    - Gets user defined data
    - Joins dataframes, replacing international data with user provided data if avail
    """

    df = costs.copy()
    df["YEAR"] = df.YEAR.astype(int)

    cmo = _get_cmo_template(df, mining_fuels, countries, years, "MIN")
    user = _get_user_template(df, mining_fuels, years, "MIN")
    final = pd.concat([user, cmo]).drop_duplicates(keep="first")

    # format data

    final["REGION"] = region
    final["MODE_OF_OPERATION"] = 1

    return final.set_index(["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"])


def get_international_mining_data(
    costs: pd.DataFrame,
    mining_fuels: list[str],
    countries: list[str],
    years: list[int],
    region: str,
    multiplier: float,
    fuel_limits: pd.DataFrame,
) -> pd.DataFrame:
    """Gets domestic mining data (ie. operates on mode 2)"""

    df = get_domestic_mining_data(
        costs, mining_fuels, countries, years, region
    ).reset_index()
    df["MODE_OF_OPERATION"] = 2
    df["VALUE"] = df.VALUE.mul(multiplier)

    df = assign_international_limits(df, fuel_limits)

    return df.set_index(["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"])


def get_domestic_renewable_data(
    costs: pd.DataFrame,
    renewable_fuels: list[str],
    nodes: list[str],
    years: list[int],
    region: str,
) -> pd.DataFrame:
    """Gets domestic mining data (ie. operates on mode 1)

    Logic is:
    - Assigns data for all regions based on international data
    - Gets user defined data
    - Joins dataframes, replacing international data with user provided data if avail
    - Explodes country data based on node (ie. all nodes within a counrty have same price)
    """

    def _get_country_nodes_mapper(nodes: list[str]) -> dict[str, list[str]]:
        country = [x[:3] for x in nodes]
        node = [x[3:] for x in nodes]

        data = {}

        for c, n in zip(country, node):
            if c in data:
                data[c].append(n)
            else:
                data[c] = [n]

        return data

    df = costs.copy()
    df["YEAR"] = df.YEAR.astype(int)
    countries_nodes = _get_country_nodes_mapper(nodes)
    countries = countries_nodes.keys()

    cmo = _get_cmo_template(df, renewable_fuels, countries, years, "RNW")
    user = _get_user_template(df, renewable_fuels, years, "RNW")
    final = pd.concat([user, cmo]).drop_duplicates(keep="first")

    # format data
    final["COUNTRY"] = final.TECHNOLOGY.str[6:]
    final["NODE"] = final.COUNTRY.map(countries_nodes)
    final = final.dropna(subset="NODE").explode("NODE")

    final.TECHNOLOGY = final.TECHNOLOGY + final.NODE

    final = final.drop(columns=["COUNTRY", "NODE"])

    final["REGION"] = region
    final["MODE_OF_OPERATION"] = 1

    return final.set_index(["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"])


def main(
    cmo_forecasts: pd.DataFrame,
    cmo_data_year: int,
    user_fuel_prices: pd.DataFrame,
    technologies: pd.Series,
    years: pd.Series,
    regions: pd.Series,
    biomass_costs: int | float,
    nuclear_costs: int | float,
    waste_costs: int | float,
    international_cost_factor: Optional[int | float] = None,
    fuel_limits: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Creates variable cost data"""

    cmo = calculate_cmo_forecasts(cmo_forecasts, cmo_data_year)

    # Add in other fuels that use same fuel as CMO
    # - Cogen is powered by gas
    # - Other petroleum products are similar to oil
    # - Petroleum products are similar to oil
    cog = copy_cmo_fuel_price(cmo, "GAS", "COG")
    oth = copy_cmo_fuel_price(cmo, "OIL", "OTH")
    pet = copy_cmo_fuel_price(cmo, "OIL", "PET")

    # Add hardcoded fuel prices
    bio = get_static_fuel_price("BIO", "INT", "m$/PJ", 1, biomass_costs)
    was = get_static_fuel_price("WAS", "INT", "m$/PJ", 1, waste_costs)
    nuc = get_static_fuel_price("URN", "INT", "m$/PJ", 1, nuclear_costs)

    # collect all cmo formatted data
    cmo = pd.concat([cmo, cog, oth, pet, bio, was, nuc])

    # assign cmo data the same years as user defined prices
    user_years = get_user_fuel_years(user_fuel_prices)
    cmo = expand_cmo_data(cmo, user_years)

    # combine cmo and user data
    df = merge_cmo_user(cmo, user_fuel_prices)

    # scaffolding information
    countries = get_countries_from_techs(technologies)
    nodes = get_nodes_from_techs(technologies)
    y = years.unique().tolist()
    r = regions.unique().tolist()[0]  # only one region

    mining_fuels = MINING_FUELS
    renewable_fuels = RENEWABLE_FUELS

    if not isinstance(fuel_limits, pd.DataFrame):
        fuel_limits = pd.DataFrame()  # no international trading

    # process cost data
    df = apply_energy_content(df)

    mining_domestic = get_domestic_mining_data(df, mining_fuels, countries, y, r)
    mining_international = get_international_mining_data(
        df, mining_fuels, countries, y, r, international_cost_factor, fuel_limits
    )
    renewable_domestic = get_domestic_renewable_data(df, renewable_fuels, nodes, y, r)

    var_costs = pd.concat([mining_domestic, mining_international, renewable_domestic])
    var_costs = filter_var_cost_technologies(var_costs, technologies)
    return var_costs.round(3)


if __name__ == "__main__":

    if "snakemake" in globals():
        file_cmo_forecasts = snakemake.input.cmo_forecasts
        file_fuel_prices = snakemake.input.fuel_prices
        file_regions = snakemake.input.regions
        file_years = snakemake.input.years
        file_technologies = snakemake.input.technologies
        file_fuel_limits = snakemake.input.fuel_limits
        file_var_costs = snakemake.output.var_costs
    else:
        file_cmo_forecasts = "resources/data/CMO-October-2024-Forecasts.xlsx"
        file_fuel_prices = "resources/data/fuel_prices.csv"
        file_regions = "results/India/data/REGION.csv"
        file_years = "results/India/data/YEAR.csv"
        file_technologies = "results/India/data/TECHNOLOGY.csv"
        file_fuel_limits = "resources/data/fuel_limits.csv"
        file_var_costs = "results/India/data/VariableCosts.csv"

    cmo_forecasts = import_cmo_forecasts(file_cmo_forecasts)
    user_fuel_prices = import_fuel_prices(file_fuel_prices)
    user_fuel_limits = import_fuel_limits(file_fuel_limits)
    technologies = import_set(file_technologies)
    years = import_set(file_years)
    regions = import_set(file_regions)

    input_data = {
        "cmo_forecasts": cmo_forecasts,
        "user_fuel_prices": user_fuel_prices,
        "technologies": technologies,
        "years": years,
        "regions": regions,
        "cmo_data_year": CMO_DATA_YEAR,
        "biomass_costs": BIOMASS_VAR_COSTS,
        "nuclear_costs": NUCLEAR_VAR_COSTS,
        "waste_costs": WASTE_VAR_COSTS,
        "international_cost_factor": INT_COST_FACTOR,
        "fuel_limits": user_fuel_limits,
    }

    df = main(**input_data)

    df.to_csv(file_var_costs, index=True)
