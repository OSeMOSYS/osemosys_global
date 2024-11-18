"""Calcualtes carbon intensity

Since emissions are reported at country levels, carbon intensities are also 
reported at country levels. 
"""

import pandas as pd


def format_demand(demand: pd.DataFrame) -> pd.DataFrame:

    df = demand.copy()

    df["CTRY"] = df.FUEL.str[3:6]
    df = df.drop(columns=["TIMESLICE", "FUEL"])
    return df.groupby(["REGION", "CTRY", "YEAR"]).sum()


def format_emissions(annual_emissions: pd.DataFrame) -> pd.DataFrame:

    df = annual_emissions.copy()

    df["CTRY"] = df.EMISSION.str[3:6]
    return df.groupby(["REGION", "EMISSION", "CTRY", "YEAR"]).sum()


def calculate_emission_intensity(
    demand: pd.DataFrame, emissions: pd.DataFrame
) -> pd.DataFrame:

    d = demand.copy().rename(columns={"VALUE": "DEMAND_PJ"})
    e = emissions.copy().rename(columns={"VALUE": "EMISSIONS_MT"})

    df = d.join(e).fillna(0)
    df = df.droplevel("CTRY")  # emission retains country

    # 1 PJ (1000 TJ / PJ) (1000 GJ / TJ) (1000 MJ / GJ) (1000 kJ / MJ) (hr / 3600sec)
    df["DEMAND_kwh"] = df.DEMAND_PJ.mul(1000).mul(1000).mul(1000).mul(1000).div(3600)

    # 1 MT (1,000,000 T / MT) (1000 kg / T) (1000g / kg)
    df["EMISSIONS_g"] = df.EMISSIONS_MT.mul(1000000).mul(1000).mul(1000)

    # intensity in g/kwh
    df["VALUE"] = df.EMISSIONS_g.div(df.DEMAND_kwh)

    return df["VALUE"].to_frame()


if __name__ == "__main__":
    if "snakemake" in globals():
        demand_csv = snakemake.input.demand
        annual_emissions_csv = snakemake.input.annual_emissions
        save = snakemake.output.emission_intensity
    else:
        demand_csv = "results/India/results/Demand.csv"
        annual_emissions_csv = "results/India/results/AnnualEmissions.csv"
        save = "results/India/results/AnnualEmissionIntensity.csv"

    demand = pd.read_csv(demand_csv)
    annual_emissions = pd.read_csv(annual_emissions_csv)

    demand = format_demand(demand)
    emissions = format_emissions(annual_emissions)

    emission_intensity = calculate_emission_intensity(demand, emissions).round(2)

    emission_intensity.to_csv(save, index=True)
