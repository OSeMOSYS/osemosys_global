"""Calcualtes carbon intensity

Since emissions are reported at country levels, carbon intensities are also 
reported at country levels. 
"""

import pandas as pd
from typing import Optional


def format_production(
    production: pd.DataFrame,
    exclude: Optional[list[str]] = None,
) -> pd.DataFrame:

    df = production.copy()

    df = df[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("PWR"))
        & ~(df.index.get_level_values("TECHNOLOGY").str.contains("TRN"))
    ].copy()

    df["TECH"] = df.index.get_level_values("TECHNOLOGY").str[3:6]
    df["COUNTRY"] = df.index.get_level_values("TECHNOLOGY").str[6:9]

    if exclude:
        df = df[~df.TECH.isin(exclude)].copy()

    return (
        df.reset_index()[["REGION", "COUNTRY", "YEAR", "VALUE"]]
        .groupby(["REGION", "COUNTRY", "YEAR"])
        .sum()
    )


def format_emissions(annual_emissions: pd.DataFrame) -> pd.DataFrame:

    df = annual_emissions.copy().reset_index()

    df["COUNTRY"] = df.EMISSION.str[3:6]
    return df.groupby(["REGION", "EMISSION", "COUNTRY", "YEAR"]).sum()


def calculate_emission_intensity(
    production: pd.DataFrame, emissions: pd.DataFrame
) -> pd.DataFrame:

    p = production.copy().rename(columns={"VALUE": "PROD_PJ"})
    e = emissions.copy().rename(columns={"VALUE": "EMISSIONS_MT"})

    df = p.join(e).fillna(0)
    df = df.droplevel("COUNTRY")  # emission retains country

    # 1 PJ (1000 TJ / PJ) (1000 GJ / TJ) (1000 MJ / GJ) (1000 kJ / MJ) (hr / 3600sec)
    df["PROD_kwh"] = df.PROD_PJ.mul(1000).mul(1000).mul(1000).mul(1000).div(3600)

    # 1 MT (1,000,000 T / MT) (1000 kg / T) (1000g / kg)
    df["EMISSIONS_g"] = df.EMISSIONS_MT.mul(1000000).mul(1000).mul(1000)

    # intensity in g/kwh
    df["VALUE"] = df.EMISSIONS_g.div(df.PROD_kwh)

    return df["VALUE"].to_frame()


if __name__ == "__main__":
    if "snakemake" in globals():
        production_csv = snakemake.input.production_by_technology
        annual_emissions_csv = snakemake.input.annual_emissions
        save = snakemake.output.emission_intensity
        storage = snakemake.params.storage
    else:
        production_csv = "results/India/results/ProductionByTechnology.csv"
        annual_emissions_csv = "results/India/results/AnnualEmissions.csv"
        save = "results/India/results/AnnualEmissionIntensity.csv"
        storage = {"SDS": [], "LDS": []}

    exclusions = list(storage)

    production = pd.read_csv(production_csv, index_col=[0, 1, 2, 3, 4])
    annual_emissions = pd.read_csv(annual_emissions_csv, index_col=[0, 1, 2])

    production = format_production(production, exclusions)
    emissions = format_emissions(annual_emissions)

    emission_intensity = calculate_emission_intensity(production, emissions).round(2)

    emission_intensity.to_csv(save, index=True)
