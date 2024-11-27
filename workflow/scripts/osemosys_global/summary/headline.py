"""Calculates Headline Metrics"""

import pandas as pd
from typing import Optional
from constants import RENEWABLES, FOSSIL, CLEAN


def get_emissions(annual_emissions: pd.DataFrame) -> pd.DataFrame:

    df = annual_emissions.copy()

    total_emissions = df.VALUE.sum().round(0)
    data = ["Emissions", "Million tonnes of CO2-eq.", total_emissions]
    return pd.DataFrame([data], columns=["Metric", "Unit", "Value"])


def get_system_cost(total_discounted_cost: pd.DataFrame) -> pd.DataFrame:

    df = total_discounted_cost.copy()

    total_cost = round((df.VALUE.sum() / 1000), 0)

    data = ["Total system cost", "Billion $", total_cost]
    return pd.DataFrame([data], columns=["Metric", "Unit", "Value"])


def get_gen_cost(
    total_discounted_cost: pd.DataFrame, demand: pd.DataFrame
) -> pd.DataFrame:

    costs = total_discounted_cost.copy()
    dem = demand.copy()

    total_cost = costs.VALUE.sum()  # $M
    total_demand = dem.VALUE.sum()  # PJ

    # ($M / PJ) ($1000000 / $M) (1 PJ / 1000 TJ) (1 TJ / 1000 GJ) (1 GJ / 1000 MJ) (3600sec / 1000 hr) = 3.6
    gen_cost = ((total_cost / total_demand) * 3.6).round(0)

    data = ["Cost of electricity", "$/MWh", gen_cost]
    return pd.DataFrame([data], columns=["Metric", "Unit", "Value"])


def _filter_pwr_techs(
    production_by_technology: pd.DataFrame, exclusions: Optional[list[str]]
) -> pd.DataFrame:
    """Filters for only power techs

    exclusions: Optional[list[str]]
        Additional 'PWR' techs to filter for. For example, if ['LDS' and 'SDS'] are provided,
        will filter out 'PWRLDS' and 'PWRSDS' values
    """

    df = production_by_technology.copy()

    df = df[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("PWR"))
        & ~(df.index.get_level_values("TECHNOLOGY").str.contains("TRN"))
    ].copy()

    if not exclusions:
        return df

    prefixes = tuple([f"PWR{x}" for x in exclusions])

    return df[~df.index.get_level_values("TECHNOLOGY").str.startswith(prefixes)].copy()


def _filter_techs(
    production_by_technology: pd.DataFrame, carriers: list[str]
) -> pd.DataFrame:

    df = production_by_technology.copy()
    df["tech"] = df.index.get_level_values("TECHNOLOGY").str[0:6]

    techs = [f"PWR{x}" for x in carriers]

    df = df[df.tech.isin(techs)]

    return df["VALUE"].to_frame()


def get_gen_shares(
    production_by_technology: pd.DataFrame, exclusions: Optional[list[str]] = None
) -> pd.DataFrame:

    gen_total = _filter_pwr_techs(production_by_technology, exclusions).VALUE.sum()
    rnw_total = _filter_techs(production_by_technology, RENEWABLES).VALUE.sum()
    fsl_total = _filter_techs(production_by_technology, FOSSIL).VALUE.sum()
    cln_total = _filter_techs(production_by_technology, CLEAN).VALUE.sum()

    rnw_share = round((rnw_total / gen_total) * 100, 0)
    fsl_share = round((fsl_total / gen_total) * 100, 0)
    cln_share = round((cln_total / gen_total) * 100, 0)

    data = [
        ["Fossil energy share", "%", fsl_share],
        ["Renewable energy share", "%", rnw_share],
        ["Clean energy share", "%", cln_share],
    ]

    return pd.DataFrame(data, columns=["Metric", "Unit", "Value"])


if __name__ == "__main__":
    if "snakemake" in globals():
        storage = snakemake.params.storage
        annual_emissions_csv = snakemake.input.annual_emissions
        production_by_technology_csv = snakemake.input.production_by_technology
        total_discounted_cost_csv = snakemake.input.total_discounted_cost
        demand_csv = snakemake.input.demand
        save = snakemake.output.metrics
    else:
        storage = {"SDS": [], "LDS": []}
        annual_emissions_csv = "results/India/results/AnnualEmissions.csv"
        production_by_technology_csv = (
            "results/India/results/ProductionByTechnologyAnnual.csv"
        )
        total_discounted_cost_csv = "results/India/results/TotalDiscountedCost.csv"
        demand_csv = "results/India/results/Demand.csv"
        save = "results/India/result_summaries/Metrics.csv"

    annual_emissions = pd.read_csv(annual_emissions_csv, index_col=[0, 1, 2])
    production_by_technology = pd.read_csv(
        production_by_technology_csv, index_col=[0, 1, 2, 3]
    )
    total_discounted_cost = pd.read_csv(total_discounted_cost_csv, index_col=[0, 1])
    demand = pd.read_csv(demand_csv, index_col=[0, 1, 2, 3])

    if storage:
        exclusions = list(storage)
    else:
        exclusions = []

    dfs = []

    dfs.append(get_emissions(annual_emissions))
    dfs.append(get_system_cost(total_discounted_cost))
    dfs.append(get_gen_cost(total_discounted_cost, demand))
    dfs.append(get_gen_shares(production_by_technology, exclusions))

    df = pd.concat(dfs)

    df.to_csv(save, index=False)
