"""Calculates Headline Metrics"""

import pandas as pd
from typing import Optional
from constants import RENEWABLES, FOSSIL, CLEAN
from otoole import read
from utils import get_discount_factor
from pathlib import Path


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


def get_gen_cost(cost: pd.DataFrame, demand: pd.DataFrame) -> pd.DataFrame:

    costs = cost.copy()
    dem = demand.copy()

    total_cost = costs.VALUE.sum()  # $M
    total_demand = dem.VALUE.sum()  # PJ

    # ($M / PJ) ($1000000 / $M) (1 PJ / 1000 TJ) (1 TJ / 1000 GJ) (1 GJ / 1000 MJ) (3600sec / 1000 hr) = 3.6
    gen_cost = ((total_cost / total_demand) * 3.6).round(0)

    data = ["Discounted Cost of electricity", "$/MWh", gen_cost]
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
    df = df.loc[df.index.get_level_values("FUEL").str.startswith("ELC")]

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
    df = df.loc[df.index.get_level_values("FUEL").str.startswith("ELC")]
    df["tech"] = df.index.get_level_values("TECHNOLOGY").str[0:6]

    techs = [f"PWR{x}" for x in carriers]

    df = df[df.tech.isin(techs)]

    return df["VALUE"].to_frame()


def get_gen_shares(
    production_by_technology: pd.DataFrame, exclusions: Optional[list[str]] = None
) -> pd.DataFrame:

    df = production_by_technology.copy()
    df = df.loc[df.index.get_level_values("FUEL").str.startswith("ELC")]

    gen_total = _filter_pwr_techs(df, exclusions).VALUE.sum()
    rnw_total = _filter_techs(df, RENEWABLES).VALUE.sum()
    fsl_total = _filter_techs(df, FOSSIL).VALUE.sum()
    cln_total = _filter_techs(df, CLEAN).VALUE.sum()

    rnw_share = round((rnw_total / gen_total) * 100, 2)
    fsl_share = round((fsl_total / gen_total) * 100, 2)
    cln_share = round((cln_total / gen_total) * 100, 2)

    data = [
        ["Fossil energy share", "%", fsl_share],
        ["Renewable energy share", "%", rnw_share],
        ["Clean energy share", "%", cln_share],
    ]

    return pd.DataFrame(data, columns=["Metric", "Unit", "Value"])


def _read_result_data(result_dir: str) -> dict[str, pd.DataFrame]:
    """Reads in result CSVs

    Issue with otoole config
    """
    data = {}
    files = [Path(x) for x in Path(result_dir).iterdir()]
    for f in files:
        df = pd.read_csv(f)
        df = df.set_index(df.columns[:-1].tolist())
        data[f.stem] = df
    return data


if __name__ == "__main__":
    if "snakemake" in globals():
        storage = snakemake.params.storage
        input_data_dir = snakemake.params.input_dir
        result_data_dir = snakemake.params.result_dir
        save = snakemake.output.metrics
        otoole_config = snakemake.params.otoole_input
    else:
        storage = {"SDS": [], "LDS": []}
        input_data_dir = "results/India/data"
        result_data_dir = "results/India/results"
        save = "results/India/result_summaries/Metrics.csv"
        otoole_config = "resources/otoole.yaml"

    input_data, input_defaults = read(otoole_config, "csv", input_data_dir)
    result_data = _read_result_data(result_data_dir)

    if input_data["DiscountRate"].empty:
        regions = input_data["REGION"].VALUE.to_list()
        assert len(regions) == 1
        input_data["DiscountRate"] = pd.DataFrame(
            [[regions[0], input_defaults["DiscountRate"]]], columns=["REGION", "VALUE"]
        ).set_index("REGION")

    if storage:
        exclusions = list(storage)
    else:
        exclusions = []

    years = input_data["YEAR"].VALUE.to_list()
    discount_factor = get_discount_factor(years, input_data["DiscountRate"])
    result_data["DiscountedDemand"] = result_data["Demand"].div(discount_factor)

    dfs = []

    dfs.append(get_emissions(result_data["AnnualEmissions"]))
    dfs.append(get_system_cost(result_data["TotalDiscountedCost"]))
    dfs.append(
        get_gen_cost(
            result_data["TotalDiscountedCost"], result_data["DiscountedDemand"]
        )
    )
    dfs.append(get_gen_shares(result_data["ProductionByTechnologyAnnual"], exclusions))

    df = pd.concat(dfs)

    df.to_csv(save, index=False)
