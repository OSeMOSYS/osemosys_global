"""Generates variable costs"""

import logging

logger = logging.getLogger(__name__)

import pandas as pd
from typing import Optional

from read import (
    import_cmo_forecasts,
    import_fuel_prices,
    import_set,
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


def expand_cmo_data(
    cmo: pd.DataFrame,
    years: list[int],
    countries: list[str],
    international_cost_multiplier: Optional[int | float] = None,
) -> pd.DataFrame:
    """Expands CMO data to match user defined years and all countries"""

    df_template = cmo.copy()

    for year in years:
        df_template[str(year)] = df_template["VALUE"]
    df_template = df_template.drop(columns="VALUE").set_index(
        ["FUEL", "COUNTRY", "UNIT", "ENERGY_CONTENT"]
    )

    df_international = df_template.copy()
    if international_cost_multiplier:
        df_international = df_international.mul(international_cost_multiplier)

    df_country = df_template.copy().droplevel("COUNTRY")
    df_country["COUNTRY"] = [countries] * len(df_country)
    df_country = df_country.explode("COUNTRY")
    df_country = df_country.reset_index().set_index(
        ["FUEL", "COUNTRY", "UNIT", "ENERGY_CONTENT"]
    )

    return pd.concat([df_international, df_country]).reset_index()


def expand_merged_data(costs: pd.DataFrame, modelled_years: list[int]) -> pd.DataFrame:
    """Expands cmo/user merged data to interpolate values over all model years"""

    df = costs.copy()
    df["YEAR"] = df.YEAR.astype(int)

    # get year bounds as interpolation can use years outside model scope
    min_year = min(min(modelled_years), min(costs.YEAR.astype(int)))
    max_year = max(max(modelled_years), max(costs.YEAR.astype(int)))

    df = df.set_index(["FUEL", "COUNTRY", "YEAR"])

    # new index parameters
    fules = df.index.get_level_values("FUEL").unique()
    years = range(min_year, max_year + 1)
    countries = df.index.get_level_values("COUNTRY").unique()

    # interpolate over years
    idx = pd.MultiIndex.from_product(
        [fules, countries, years], names=["FUEL", "COUNTRY", "YEAR"]
    )
    df = df.reindex(idx)
    df = df.groupby(level=["FUEL", "COUNTRY"])["VALUE"].apply(
        lambda x: x.interpolate(method="linear", limit_direction="both")
    )

    # drop level to get rid of the grouping key
    return df.droplevel([0, 1]).reset_index()


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


def get_mining_data(
    costs: pd.DataFrame, mining_fuels: list[str], region: str
) -> pd.DataFrame:
    """Gets mining data (ie. operates on mode 1)"""

    df = costs.copy()

    df = df[df.FUEL.isin(mining_fuels)].copy()

    df["TECHNOLOGY"] = "MIN" + df.FUEL + df.COUNTRY

    # format data

    df["REGION"] = region
    df["MODE_OF_OPERATION"] = 1

    return df[["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR", "VALUE"]].set_index(
        ["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"]
    )


def get_renewable_data(
    costs: pd.DataFrame,
    renewable_fuels: list[str],
    nodes: list[str],
    region: str,
) -> pd.DataFrame:
    """Gets renewable data (ie. operates on mode 1)"""

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

    df = df[df.FUEL.isin(renewable_fuels)].copy()

    countries_nodes = _get_country_nodes_mapper(nodes)

    # format data
    df["NODE"] = df.COUNTRY.map(countries_nodes)
    df = df.dropna(subset="NODE").explode("NODE")

    # situation with no user defined renewable var costs
    if df.empty:
        return pd.DataFrame(
            columns=["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR", "VALUE"]
        ).set_index(["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"])

    df["TECHNOLOGY"] = "RNW" + df.FUEL + df.COUNTRY + df.NODE

    df = df.drop(columns=["COUNTRY", "NODE", "FUEL"])

    df["REGION"] = region
    df["MODE_OF_OPERATION"] = 1

    return df.set_index(["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"])


def get_backstop_var_costs(
    techs: pd.Series, years: pd.Series, regions: pd.Series
) -> pd.Series:
    """Gets backstop variable costs"""

    t = techs[techs.str.startswith("PWRBCK")].unique().tolist()
    y = years.unique().tolist()
    r = regions.unique().tolist()[0]  # only one region

    df = pd.DataFrame(
        index=pd.MultiIndex.from_product(
            [[r], t, y], names=["REGION", "TECHNOLOGY", "YEAR"]
        )
    ).reset_index()
    df["MODE_OF_OPERATION"] = 1
    df["VALUE"] = 999999

    return df.set_index(["REGION", "TECHNOLOGY", "MODE_OF_OPERATION", "YEAR"])


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
) -> pd.DataFrame:
    """Creates variable cost data"""

    # scaffolding information
    countries = get_countries_from_techs(technologies)
    nodes = get_nodes_from_techs(technologies)
    y = years.unique().tolist()
    r = regions.unique().tolist()[0]  # only one region

    # international markets data
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
    # note cost multipler must be applied here, as we do not want to modify
    # the user defined data with the multiplier
    user_years = get_user_fuel_years(user_fuel_prices)
    cmo = expand_cmo_data(cmo, user_years, countries, international_cost_factor)

    # combine cmo and user data
    # user data is used wherever possible, with cmo to fill in gaps
    df = merge_cmo_user(cmo, user_fuel_prices)

    mining_fuels = MINING_FUELS
    renewable_fuels = RENEWABLE_FUELS

    # process cost data
    df = apply_energy_content(df)
    df = expand_merged_data(df, years)

    mining_data = get_mining_data(df, mining_fuels, r)
    renewable_data = get_renewable_data(df, renewable_fuels, nodes, r)

    var_costs = pd.concat([mining_data, renewable_data])
    var_costs = var_costs[var_costs.index.get_level_values("YEAR").isin(y)]
    var_costs = filter_var_cost_technologies(var_costs, technologies)

    bck_var_costs = get_backstop_var_costs(technologies, years, regions)

    var_costs = pd.concat([var_costs, bck_var_costs])

    return var_costs.round(3)


if __name__ == "__main__":

    if "snakemake" in globals():
        file_cmo_forecasts = snakemake.input.cmo_forecasts
        file_fuel_prices = snakemake.input.fuel_prices
        file_regions = snakemake.input.regions
        file_years = snakemake.input.years
        file_technologies = snakemake.input.technologies
        file_var_costs = snakemake.output.var_costs
    else:
        file_cmo_forecasts = "resources/data/CMO-October-2024-Forecasts.xlsx"
        file_fuel_prices = "resources/data/fuel_prices.csv"
        file_regions = "results/data/REGION.csv"
        file_years = "results/data/YEAR.csv"
        file_technologies = "results/data/powerplant/TECHNOLOGY.csv"
        file_var_costs = "results/data/powerplant/VariableCosts.csv"

    cmo_forecasts = import_cmo_forecasts(file_cmo_forecasts)
    user_fuel_prices = import_fuel_prices(file_fuel_prices)
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
    }

    df = main(**input_data)

    df.to_csv(file_var_costs, index=True)
