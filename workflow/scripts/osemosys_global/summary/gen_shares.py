"""Calcualtes Renewable Generation"""

import pandas as pd
from typing import Optional
from constants import CLEAN, RENEWABLES, FOSSIL


def _get_gen_by_node(
    production: pd.DataFrame,
    include: Optional[list[str]] = None,
    exclude: Optional[list[str]] = None,
) -> pd.DataFrame:

    df = production.copy()
    df = df.loc[df.index.get_level_values("FUEL").str.startswith('ELC')]

    assert "TECHNOLOGY" in df.index.names

    df["TECH"] = df.index.get_level_values("TECHNOLOGY").str[3:6]
    df["NODE"] = df.index.get_level_values("TECHNOLOGY").str[6:11]

    if exclude:
        df = df[~df.TECH.isin(exclude)].copy()

    if include:
        df = df[df.TECH.isin(include)].copy()

    return (
        df.reset_index()[["REGION", "NODE", "YEAR", "VALUE"]]
        .groupby(["REGION", "NODE", "YEAR"])
        .sum()
    )


def _get_gen_by_country(
    production: pd.DataFrame,
    include: Optional[list[str]] = None,
    exclude: Optional[list[str]] = None,
) -> pd.DataFrame:

    df = production.copy()
    df = df.loc[df.index.get_level_values("FUEL").str.startswith('ELC')]

    assert "TECHNOLOGY" in df.index.names

    df["TECH"] = df.index.get_level_values("TECHNOLOGY").str[3:6]
    df["COUNTRY"] = df.index.get_level_values("TECHNOLOGY").str[6:9]

    if exclude:
        df = df[~df.TECH.isin(exclude)].copy()

    if include:
        df = df[df.TECH.isin(include)].copy()

    return (
        df.reset_index()[["REGION", "COUNTRY", "YEAR", "VALUE"]]
        .groupby(["REGION", "COUNTRY", "YEAR"])
        .sum()
    )

def _get_gen_global(
    production: pd.DataFrame,
    include: Optional[list[str]] = None,
    exclude: Optional[list[str]] = None,
) -> pd.DataFrame:

    df = production.copy()
    df = df.loc[df.index.get_level_values("FUEL").str.startswith('ELC')]

    assert "TECHNOLOGY" in df.index.names

    df["TECH"] = df.index.get_level_values("TECHNOLOGY").str[3:6]
    #df["COUNTRY"] = df.index.get_level_values("TECHNOLOGY").str[6:9]

    if exclude:
        df = df[~df.TECH.isin(exclude)].copy()

    if include:
        df = df[df.TECH.isin(include)].copy()

    return (
        df.reset_index()[["REGION", "YEAR", "VALUE"]]
        .groupby(["REGION", "YEAR"])
        .sum()
    )

def calc_generation_shares_node(
    production_by_technology: pd.DataFrame, exclusions: Optional[list[str]] = None
) -> pd.DataFrame:

    if not exclusions:
        exclusions = []

    df = production_by_technology.copy()

    df = df[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("PWR"))
        & ~(df.index.get_level_values("TECHNOLOGY").str.contains("TRN"))
    ]

    total = _get_gen_by_node(df, exclude=exclusions).rename(columns={"VALUE": "TOTAL"})
    
    clean = _get_gen_by_node(df, include=CLEAN).rename(columns={"VALUE": "CLEAN"})
    
    renewable = _get_gen_by_node(df, include=RENEWABLES).rename(
        columns={"VALUE": "RENEWABLE"}
    )
    
    fossil = _get_gen_by_node(df, include=FOSSIL).rename(columns={"VALUE": "FOSSIL"})

    shares = (
        total.join(clean, how="outer")
        .join(renewable, how="outer")
        .join(fossil, how="outer")
        .fillna(0)
    )

    shares["CLEAN"] = shares.CLEAN.div(shares.TOTAL).mul(100)
    shares["RENEWABLE"] = shares.RENEWABLE.div(shares.TOTAL).mul(100)
    shares["FOSSIL"] = shares.FOSSIL.div(shares.TOTAL).mul(100)

    return shares[["CLEAN", "RENEWABLE", "FOSSIL"]].round(1)


def calc_generation_shares_country(
    production_by_technology: pd.DataFrame, exclusions: Optional[list[str]] = None
) -> pd.DataFrame:

    if not exclusions:
        exclusions = []

    df = production_by_technology.copy()

    df = df[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("PWR"))
        & ~(df.index.get_level_values("TECHNOLOGY").str.contains("TRN"))
    ]

    total = _get_gen_by_country(df, exclude=exclusions).rename(
        columns={"VALUE": "TOTAL"}
    )
    clean = _get_gen_by_country(df, include=CLEAN).rename(columns={"VALUE": "CLEAN"})
    renewable = _get_gen_by_country(df, include=RENEWABLES).rename(
        columns={"VALUE": "RENEWABLE"}
    )
    fossil = _get_gen_by_country(df, include=FOSSIL).rename(columns={"VALUE": "FOSSIL"})

    shares = (
        total.join(clean, how="outer")
        .join(renewable, how="outer")
        .join(fossil, how="outer")
        .fillna(0)
    )

    shares["CLEAN"] = shares.CLEAN.div(shares.TOTAL).mul(100)
    shares["RENEWABLE"] = shares.RENEWABLE.div(shares.TOTAL).mul(100)
    shares["FOSSIL"] = shares.FOSSIL.div(shares.TOTAL).mul(100)

    return shares[["CLEAN", "RENEWABLE", "FOSSIL"]].round(1)

def calc_generation_shares_global(
    production_by_technology: pd.DataFrame, exclusions: Optional[list[str]] = None
) -> pd.DataFrame:

    if not exclusions:
        exclusions = []

    df = production_by_technology.copy()

    df = df[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("PWR"))
        & ~(df.index.get_level_values("TECHNOLOGY").str.contains("TRN"))
    ]

    total = _get_gen_global(df, exclude=exclusions).rename(
        columns={"VALUE": "TOTAL"}
    )
    clean = _get_gen_global(df, include=CLEAN).rename(columns={"VALUE": "CLEAN"})
    renewable = _get_gen_global(df, include=RENEWABLES).rename(
        columns={"VALUE": "RENEWABLE"}
    )
    fossil = _get_gen_global(df, include=FOSSIL).rename(columns={"VALUE": "FOSSIL"})

    shares = (
        total.join(clean, how="outer")
        .join(renewable, how="outer")
        .join(fossil, how="outer")
        .fillna(0)
    )

    shares["CLEAN"] = shares.CLEAN.div(shares.TOTAL).mul(100)
    shares["RENEWABLE"] = shares.RENEWABLE.div(shares.TOTAL).mul(100)
    shares["FOSSIL"] = shares.FOSSIL.div(shares.TOTAL).mul(100)

    return shares[["CLEAN", "RENEWABLE", "FOSSIL"]].round(1)


if __name__ == "__main__":
    if "snakemake" in globals():
        production_by_technology_annual_csv = snakemake.input.production_by_technology
        storage = snakemake.params.storage
        gen_shares_node = snakemake.output.generation_shares_node
        gen_shares_country = snakemake.output.generation_shares_country
        gen_shares_global = snakemake.output.generation_shares_global
    else:
        production_by_technology_annual_csv = (
            "results/India/results/ProductionByTechnologyAnnual.csv"
        )
        storage = {"SDS": [], "LDS": []}
        gen_shares_node = "results/India/result_summaries/GenerationSharesNode.csv"
        gen_shares_country = (
            "results/India/result_summaries/GenerationSharesCountry.csv"
        )
        gen_shares_global = (
            "results/India/result_summaries/GenerationSharesGlobal.csv"
        )


    production_by_technology_annual = pd.read_csv(
        production_by_technology_annual_csv, index_col=[0, 1, 2, 3]
    )

    if storage:
        exclusions = list(storage)
    else:
        exclusions = []

    nodes = calc_generation_shares_node(production_by_technology_annual, exclusions)
    country = calc_generation_shares_country(production_by_technology_annual, exclusions)
    glob = calc_generation_shares_global(production_by_technology_annual, exclusions)

    nodes.to_csv(gen_shares_node, index=True)
    country.to_csv(gen_shares_country, index=True)
    glob.to_csv(gen_shares_global, index=True)
