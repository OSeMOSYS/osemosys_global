"""Summarizes tranmsission and power capacity"""

import pandas as pd

def calc_trn_capacity(total_capacity_annual: pd.DataFrame) -> pd.DataFrame:

    df = total_capacity_annual.copy()

    df = df[df.index.get_level_values("TECHNOLOGY").str.startswith("TRN")]
    df["FROM"] = df.index.get_level_values("TECHNOLOGY").str[3:8]
    df["TO"] = df.index.get_level_values("TECHNOLOGY").str[8:13]

    return df.reset_index()[["FROM", "TO", "YEAR", "VALUE"]].groupby(["FROM", "TO", "YEAR"]).sum()

def calc_pwr_capacity(total_capacity_annual: pd.DataFrame) -> pd.DataFrame:

    df = total_capacity_annual.copy()

    df = df[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("PWR"))
        & ~(df.index.get_level_values("TECHNOLOGY").str.contains("TRN"))
    ]
    
    df["TECH"] = df.index.get_level_values("TECHNOLOGY").str[3:6]
    df["NODE"] = df.index.get_level_values("TECHNOLOGY").str[6:11]

    return df.reset_index()[["TECH", "NODE", "YEAR", "VALUE"]].groupby(["TECH", "NODE", "YEAR"]).sum().rename(columns={"TECH": "TECHNOLOGY"})

if __name__ == "__main__":
    if "snakemake" in globals():
        total_capacity_csv = snakemake.input.total_capacity
        pwr_save = snakemake.output.power_capacity
        trn_save = snakemake.output.transmission_capacity
    else:
        total_capacity_csv = (
            "results/India/results/TotalCapacityAnnual.csv"
        )
        pwr_save = ""
        trn_save = ""

    total_capacity_annual = pd.read_csv(
        total_capacity_csv, index_col=[0, 1, 2]
    )

    pwr_capacity = calc_pwr_capacity(total_capacity_annual)
    trn_capacity = calc_trn_capacity(total_capacity_annual)

    pwr_capacity.to_csv(pwr_save, index=True)
    trn_capacity.to_csv(trn_save, index=True)
