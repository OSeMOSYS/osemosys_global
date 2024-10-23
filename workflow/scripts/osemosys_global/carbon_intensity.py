"""Calcualtes carbon intensity

Since emissions are reported at country levels, carbon intensities are also 
reported at country levels. 
"""

import pandas as pd


def format_production(production_by_technology_annual: pd.DataFrame) -> pd.DataFrame:

    df = production_by_technology_annual.copy()

    df = df[df.TECHNOLOGY.str.startswith("PWR")]
    df["TECH"] = df.TECHNOLOGY.str[3:6]
    df["CTRY"] = df.TECHNOLOGY.str[6:9]
    df = df[~(df.TECH.isin(["BAT", "TRN"]))]
    df = df.drop(columns=["TECHNOLOGY", "FUEL", "TECH"])
    return df.groupby(["REGION", "CTRY", "YEAR"]).sum()


def format_emissions(annual_emissions: pd.DataFrame) -> pd.DataFrame:

    df = annual_emissions.copy()

    df["CTRY"] = df.EMISSION.str[3:6]
    return df.groupby(["REGION", "EMISSION", "CTRY", "YEAR"]).sum()


def calculate_emission_intensity(
    production: pd.DataFrame, emissions: pd.DataFrame
) -> pd.DataFrame:

    p = production.copy().rename(columns={"VALUE": "PRODUCTION_PJ"})
    e = emissions.copy().rename(columns={"VALUE": "EMISSIONS_MT"})

    df = p.join(e).fillna(0)
    df = df.droplevel("CTRY")  # emission retains country

    # 1 PJ (1000 TJ / PJ) (1000 GJ / TJ) (1000 MJ / GJ) (1000 kJ / MJ) (hr / 3600sec)
    df["PRODUCTION_kwh"] = (
        df.PRODUCTION_PJ.mul(1000).mul(1000).mul(1000).mul(1000).div(3600)
    )

    # 1 MT (1,000,000 T / MT) (1000 kg / T) (1000g / kg)
    df["EMISSIONS_g"] = df.EMISSIONS_MT.mul(1000000).mul(1000).mul(1000)

    # intensity in g/kwh
    df["VALUE"] = df.EMISSIONS_g.div(df.PRODUCTION_kwh)

    return df["VALUE"].to_frame()


if __name__ == "__main__":
    if "snakemake" in globals():
        production_by_technology_annual_csv = snakemake.input.production_by_technology
        annual_emissions_csv = snakemake.input.annual_emissions
        save = snakemake.output.emission_intensity
    else:
        production_by_technology_annual_csv = (
            "results/India/results/ProductionByTechnologyAnnual.csv"
        )
        annual_emissions_csv = "results/India/results/AnnualEmissions.csv"
        save = "results/India/results/AnnualEmissionIntensity.csv"

    production_by_technology_annual = pd.read_csv(production_by_technology_annual_csv)
    annual_emissions = pd.read_csv(annual_emissions_csv)

    prod = format_production(production_by_technology_annual)
    emissions = format_emissions(annual_emissions)

    emission_intensity = calculate_emission_intensity(prod, emissions).round(2)

    emission_intensity.to_csv(save, index=True)
