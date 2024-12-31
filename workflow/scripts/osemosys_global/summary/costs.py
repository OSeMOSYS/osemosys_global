"""Gets power generation cost per node"""

import pandas as pd
from pathlib import Path


def get_tech_cost(discounted_cost_tech: pd.DataFrame, country: bool) -> pd.DataFrame:
    """Only of power generation technologies"""

    df = discounted_cost_tech.copy()

    df = df[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("PWR"))
        & ~(df.index.get_level_values("TECHNOLOGY").str.contains("TRN"))
    ]

    if country:
        r = "COUNTRY"
        df[r] = df.index.get_level_values("TECHNOLOGY").str[6:9]
    else:
        r = "NODE"
        df[r] = df.index.get_level_values("TECHNOLOGY").str[6:11]

    return (
        df.reset_index()[["REGION", r, "YEAR", "VALUE"]]
        .groupby(["REGION", r, "YEAR"])
        .sum()
    )


def get_storage_cost(
    discounted_cost_storage: pd.DataFrame, country: bool
) -> pd.DataFrame:

    df = discounted_cost_storage.copy()

    if country:
        r = "COUNTRY"
        df[r] = df.index.get_level_values("STORAGE").str[3:6]
    else:
        r = "NODE"
        df[r] = df.index.get_level_values("STORAGE").str[3:8]

    return (
        df.reset_index()[["REGION", r, "YEAR", "VALUE"]]
        .groupby(["REGION", r, "YEAR"])
        .sum()
    )

def get_transmission_cost(discounted_cost_tech: pd.DataFrame, country: bool) -> pd.DataFrame:

    df = discounted_cost_tech.copy()

    df1 = df.copy()[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("TRN"))
    ]

    if country:
        r = "COUNTRY"
        df1[r] = df1.index.get_level_values("TECHNOLOGY").str[3:6]
    else:
        r = "NODE"
        df1[r] = df1.index.get_level_values("TECHNOLOGY").str[3:8]
    
    df2 = df.copy()[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("TRN"))
    ]

    if country:
        r = "COUNTRY"
        df2[r] = df2.index.get_level_values("TECHNOLOGY").str[8:11]
    else:
        r = "NODE"
        df2[r] = df2.index.get_level_values("TECHNOLOGY").str[8:13]
        
    for df in [df1, df2]:
        df['VALUE'] = df['VALUE'] / 2
        
    df = pd.concat([df1, df2])

    return (
        df.reset_index()[["REGION", r, "YEAR", "VALUE"]]
        .groupby(["REGION", r, "YEAR"])
        .sum()
    )

def get_demand(demand: pd.DataFrame, country: bool) -> pd.DataFrame:

    df = demand.copy()

    df = df[df.index.get_level_values("FUEL").str.startswith("ELC")]

    if country:
        r = "COUNTRY"
        df[r] = df.index.get_level_values("FUEL").str[3:6]
    else:
        r = "NODE"
        df[r] = df.index.get_level_values("FUEL").str[3:8]

    return (
        df.reset_index()[["REGION", r, "YEAR", "VALUE"]]
        .groupby(["REGION", r, "YEAR"])
        .sum()
    )


def get_pwr_cost(demand: pd.DataFrame, cost: pd.DataFrame) -> pd.DataFrame:
    """Gets power generation cost per node in $/MWh"""

    # ($M / PJ) (1PJ / 1000 TJ) (1TJ / 1000 GJ) (1GJ / 1000 MJ) ($1000000 / $M) (3600sec / hr)
    df = cost.div(demand)
    return df.mul(3.6)


if __name__ == "__main__":
    if "snakemake" in globals():
        discounted_cost_by_technology_csv = (
            snakemake.input.discounted_cost_by_technology
        )
        demand_csv = snakemake.input.demand
        power_cost_node_csv = snakemake.output.node_pwr_cost
        power_cost_country_csv = snakemake.output.country_pwr_cost
        power_cost_global_csv = snakemake.output.global_pwr_cost        
        total_cost_node_csv = snakemake.output.node_cost
        total_cost_country_csv = snakemake.output.country_cost
        total_cost_global_csv = snakemake.output.global_cost
        
    else:
        discounted_cost_by_technology_csv = (
            "results/India/results/DiscountedCostByTechnology.csv"
        )
        demand_csv = "results/India/results/Demand.csv"
        power_cost_node_csv = "results/India/result_summaries/PowerCostNode.csv"
        power_cost_country_csv = "results/India/result_summaries/PowerCostCountry.csv"
        power_cost_global_csv = "results/India/result_summaries/PowerCostGlobal.csv"
        total_cost_node_csv = "results/India/result_summaries/TotalCostNode.csv"
        total_cost_country_csv = "results/India/result_summaries/TotalCostCountry.csv"
        total_cost_global_csv = "results/India/result_summaries/TotalCostGlobal.csv"

    discounted_cost_by_technology = pd.read_csv(
        discounted_cost_by_technology_csv, index_col=[0, 1, 2]
    )
    demand_raw = pd.read_csv(demand_csv, index_col=[0, 1, 2, 3])

    # handle the storage discounted costs outside of snakemake
    # if running a model without storage, its not guaranteed that result
    # csvs will be generated with all solvers
    try:
        parent = Path(discounted_cost_by_technology_csv).parent
        discounted_cost_by_storage_csv = Path(parent, "DiscountedCostByStorage.csv")
        discounted_cost_by_storage = pd.read_csv(
            discounted_cost_by_storage_csv, index_col=[0, 1, 2]
        )
    except FileNotFoundError:
        discounted_cost_by_storage = pd.DataFrame(
            columns=["REGION", "STORAGE", "YEAR", "VALUE"]
        ).set_index(["REGION", "STORAGE", "YEAR"])

    # node level metrics

    tech_cost = get_tech_cost(discounted_cost_by_technology, country=False)
    storage_cost = get_storage_cost(discounted_cost_by_storage, country=False)
    transmission_cost = get_transmission_cost(discounted_cost_by_technology, country=False)
    
    cost = tech_cost.add(storage_cost, fill_value = 0
                         ).add(transmission_cost, fill_value = 0)
    
    demand = get_demand(demand_raw, country=False)
    pwr_cost = get_pwr_cost(demand, cost)

    pwr_cost.to_csv(power_cost_node_csv, index=True)
    cost.to_csv(total_cost_node_csv, index=True)

    # country level metrics

    tech_cost_country = get_tech_cost(discounted_cost_by_technology, country=True)
    storage_cost_country = get_storage_cost(discounted_cost_by_storage, country=True)
    transmission_cost_country = get_transmission_cost(discounted_cost_by_technology, country=True)    
    
    cost_country = tech_cost_country.add(storage_cost_country, fill_value = 0
                                         ).add(transmission_cost_country, fill_value = 0)
    
    demand_country = get_demand(demand_raw, country=True)
    pwr_cost_country = get_pwr_cost(demand_country, cost_country)

    pwr_cost_country.to_csv(power_cost_country_csv, index=True)
    cost_country.to_csv(total_cost_country_csv, index=True)
    
    # global level metrics
    
    cost_global = tech_cost_country.add(storage_cost_country, 
                                        fill_value = 0
                                        ).add(transmission_cost_country, fill_value = 0
                                              ).groupby(["YEAR"]).sum()
  
    demand_global = get_demand(demand_raw, country=True)
    demand_global = demand_global.groupby([
        demand_global.index.get_level_values("YEAR")]).sum()

    pwr_cost_global = get_pwr_cost(demand_global, cost_global)

    pwr_cost_global.to_csv(power_cost_global_csv, index=True)
    cost_global.to_csv(total_cost_global_csv, index=True)