"""Gets power generation cost per node"""

import pandas as pd

def get_cost_per_node(discounted_cost_tech: pd.DataFrame) -> pd.DataFrame:
    """Only of power generation technologies"""

    df = discounted_cost_tech.copy()
    
    df = df[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("PWR"))
        & ~(df.index.get_level_values("TECHNOLOGY").str.contains("TRN"))
    ]
    
    df["NODE"] = df.index.get_level_values("TECHNOLOGY").str[6:11]
    
    return df.reset_index()[["REGION", "NODE", "YEAR", "VALUE"]].groupby(["REGION", "NODE", "YEAR"]).sum()
    
def get_demand_per_node(demand: pd.DataFrame) -> pd.DataFrame:
    
    df = demand.copy()
    
    df = df[df.index.get_level_values("FUEL").str.startswith("ELC")]
    df["NODE"] = df.index.get_level_values("FUEL").str[3:8]
    
    return df.reset_index()[["REGION", "NODE", "YEAR", "VALUE"]].groupby(["REGION", "NODE", "YEAR"]).sum()
    
    
def get_pwr_cost_per_node(demand: pd.DataFrame, cost: pd.DataFrame) -> pd.DataFrame:
    """Gets power generation cost per node in $/MWh"""
    
    # ($M / PJ) (1PJ / 1000 TJ) (1TJ / 1000 GJ) (1GJ / 1000 MJ) ($1000000 / $M) (3600sec / hr)
    df = cost.div(demand)
    return df.mul(3.6)
    
    
if __name__ == "__main__":
    if "snakemake" in globals():
        discounted_cost_by_technology_csv = snakemake.input.discounted_cost_by_technology
        demand_csv = snakemake.input.demand
        save = snakemake.output.node_cost
    else:
        discounted_cost_by_technology_csv = (
            "results/India/results/DiscountedCostByTechnology.csv"
        )
        demand_csv = "results/India/results/Demand.csv"
        save = "results/India/results/NodeCost.csv"

    discounted_cost_by_technology_raw = pd.read_csv(
        discounted_cost_by_technology_csv, index_col=[0, 1, 2]
    )
    demand_raw = pd.read_csv(
        demand_csv, index_col=[0, 1, 2, 3]
    )

    cost = get_cost_per_node(discounted_cost_by_technology_raw)
    demand = get_demand_per_node(demand_raw)

    df = get_pwr_cost_per_node(demand, cost)

    df.to_csv(save, index=True)
