"""Creates template generators"""

import pandas as pd
from datetime import datetime
from typing import Optional

from read import import_plexos_2015
from data import get_nodes_from_yaml, get_techs_from_yaml, Technology, Node

import logging

logger = logging.getLogger(__name__)


PLEXOS_TECH_MAPPER = {
    "Bio": "BIO",
    "Coa": "COA",
    "Cog": "COG",
    "Gas-CCGT": "CCG",
    "Gas-OCGT": "OCG",
    "Geo": "GEO",
    "Hyd": "HYD",
    "Nuc": "URN",
    "Oil": "OIL",
    "Oth": "OTH",
    "Pet": "PET",
    "Sol": "SPV",
    "Was": "WAS",
    "Wav": "WAV",
    "Win": "WON",
    "Csp": "CSP",
    "Wof": "WOF",
    "oil": "OIL",
    "gas": "OCG",
}


def get_op_life(
    techs: list[Technology], default: Optional[float] = 30
) -> dict[str, float]:
    data = {}
    for tech in techs:
        lifetime = tech.lifetime
        if not lifetime:
            data[tech.code] = default
        else:
            data[tech.code] = tech.lifetime
    return data


def get_generator_table(
    plexos_prop: pd.DataFrame,
    plexos_memb: pd.DataFrame,
    start_year: int,
    op_life_mapper: dict[str, float],
) -> pd.DataFrame:
    """Sets the main generator table derived from the PLEXOS-World model."""

    # Create main generator table
    gen_cols_1 = ["child_class", "child_object", "property", "value"]
    df_gen = plexos_prop.copy()[gen_cols_1]

    df_gen = df_gen[df_gen["child_class"] == "Generator"]
    df_gen.rename(columns={"child_object": "powerplant"}, inplace=True)
    df_gen.drop("child_class", axis=1, inplace=True)
    df_gen = pd.pivot_table(
        df_gen,
        index="powerplant",
        columns="property",
        values="value",
        aggfunc="sum",
        fill_value=0,
    )
    df_gen["total_capacity"] = (df_gen["Max Capacity"].astype(float)) * (
        df_gen["Units"].astype(int)
    )

    gen_cols_base = ["Commission Date", "Heat Rate", "Max Capacity", "total_capacity"]
    df_gen_base = df_gen[gen_cols_base]

    df_dict = plexos_memb[plexos_memb["parent_class"] == "Generator"].rename(
        {"parent_object": "powerplant"}, axis=1
    )

    ## Compile dataframe with powerplants, nodes, and fuels
    df_dict_fuel = df_dict[df_dict["collection"] == "Fuels"]
    df_dict_fuel = df_dict_fuel[["powerplant", "child_object"]]
    df_dict_nodes = df_dict[df_dict["collection"] == "Nodes"]
    df_dict_nodes = df_dict_nodes[["powerplant", "child_object"]]
    df_dict_2 = pd.merge(df_dict_fuel, df_dict_nodes, how="outer", on="powerplant")

    ## Merge original generator dataframe with nodes and fuels
    df_gen_base = pd.merge(df_gen_base, df_dict_2, how="outer", on="powerplant")
    df_gen_base.rename(
        {"child_object_x": "fuel", "child_object_y": "node"}, axis=1, inplace=True
    )

    ## Extract start year from Commission Date
    df_gen_base["Commission Date"] = pd.to_timedelta(
        df_gen_base["Commission Date"].astype(int), unit="d"
    ) + datetime(1900, 1, 1)

    df_gen_base["start_year"] = df_gen_base["Commission Date"].dt.year
    df_gen_base.drop("Commission Date", axis=1, inplace=True)

    ## Calculate efficiency from heat rate. Units of heat rate in MJ/kWh
    df_gen_base["efficiency"] = 3.6 / df_gen_base["Heat Rate"].astype(float)
    df_gen_base.drop("Heat Rate", axis=1, inplace=True)

    ## Calcluate years of operation from start year
    df_gen_base["years_of_operation"] = start_year - df_gen_base["start_year"]

    ## Fix blank spaces in 'fuels' columns. Appearing for 'Oil' powerplants in certain countries
    df_gen_base.loc[df_gen_base["fuel"].isna(), "fuel"] = (
        df_gen_base["node"].str.split("-").str[:2].str.join("-")
        + " "
        + df_gen_base["powerplant"].str.split("_", expand=True)[1]
    )

    ## Create column for technology
    df_gen_base["technology"] = df_gen_base["powerplant"].str.split("_").str[1]
    df_gen_base["technology"] = df_gen_base["technology"].str.title()

    ## Divide Gas into CCGT and OCGT based on max capacity
    df_gen_base.loc[
        (df_gen_base["technology"] == "Gas")
        & (df_gen_base["Max Capacity"].astype(float) > 130),
        "technology",
    ] = "Gas-CCGT"
    df_gen_base.loc[
        (df_gen_base["technology"] == "Gas")
        & (df_gen_base["Max Capacity"].astype(float) <= 130),
        "technology",
    ] = "Gas-OCGT"

    # Create table with aggregated capacity
    df_gen_agg_node = df_gen_base[df_gen_base["start_year"] <= start_year]
    df_gen_agg_node = df_gen_agg_node.groupby(["node", "technology"], as_index=False)[
        "total_capacity"
    ].sum()
    df_gen_agg_node = (
        df_gen_agg_node.pivot(
            index="node", columns="technology", values="total_capacity"
        )
        .fillna(0)
        .reset_index()
    )

    df_gen_agg_node.drop(
        "Sto", axis=1, inplace=True
    )  # Drop 'Sto' technology. Only for USA.

    df_gen_agg_node = (
        df_gen_agg_node.fillna(0).sort_values(by="node").set_index("node").round(2)
    )

    # Add region and country code columns
    df_gen_base["region_code"] = df_gen_base["node"].str[:2]
    df_gen_base["country_code"] = df_gen_base["node"].str[3:]

    df_gen_base["tech_code"] = df_gen_base["technology"].map(PLEXOS_TECH_MAPPER)
    df_gen_base["operational_life"] = (
        df_gen_base["tech_code"].str.lower().map(op_life_mapper)
    )

    df_gen_base["retirement_year_data"] = (
        df_gen_base["operational_life"] + df_gen_base["start_year"]
    )
    df_gen_base["retirement_diff"] = (
        df_gen_base["years_of_operation"] - df_gen_base["operational_life"]
    ) / df_gen_base["operational_life"]

    """ Set retirement year based on years of operation. 
    If (years of operation - operational life) is more than 50% of 
    operational life, set retirement year
    """
    df_gen_base.loc[df_gen_base["retirement_diff"] >= 0.5, "retirement_year_model"] = (
        2028
    )
    df_gen_base.loc[
        (df_gen_base["retirement_diff"] < 0.5) & (df_gen_base["retirement_diff"] > 0),
        "retirement_year_model",
    ] = 2033
    df_gen_base.loc[df_gen_base["retirement_diff"] <= 0, "retirement_year_model"] = (
        df_gen_base["retirement_year_data"]
    )

    df_gen_base.loc[df_gen_base["node"].str.len() <= 6, "node_code"] = (
        df_gen_base["node"].str.split("-").str[1:].str.join("") + "XX"
    )
    df_gen_base.loc[df_gen_base["node"].str.len() > 6, "node_code"] = (
        df_gen_base["node"].str.split("-").str[1:].str.join("")
    )

    df_gen_base = df_gen_base.loc[~df_gen_base["tech_code"].isna()]

    return df_gen_base


if __name__ == "__main__":
    if "snakemake" in globals():
        start_year = snakemake.params.start_year
        file_plexos = snakemake.input.plexos
        nodes_yaml = snakemake.input.nodes
        techs_yaml = snakemake.input.techs
        csv = snakemake.output
    else:
        start_year = 2021
        file_plexos = "resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx"
        nodes_yaml = "resources/data/nodes.yaml"
        techs_yaml = "resources/data/techs.yaml"
        csv = "gen_table.csv"

    plexos_prop = import_plexos_2015(file_plexos, "prop")
    plexos_memb = import_plexos_2015(file_plexos, "memb")

    techs = get_techs_from_yaml(techs_yaml)
    nodes = get_nodes_from_yaml(nodes_yaml)

    operational_life_mapper = get_op_life(techs)

    gens = get_generator_table(
        plexos_prop=plexos_prop,
        plexos_memb=plexos_memb,
        start_year=start_year,
        op_life_mapper=operational_life_mapper,
    )
    assert not gens.isnull().values.any(), "Generators contain NaN values"
    gens.to_csv(csv)
