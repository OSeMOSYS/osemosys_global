"""Summarizes tranmsission and power capacity"""

import pandas as pd


def calc_trn_capacity(
    total_capacity_annual: pd.DataFrame, country: bool
) -> pd.DataFrame:

    df = total_capacity_annual.copy()

    df = df[df.index.get_level_values("TECHNOLOGY").str.startswith("TRN")]

    if country:
        df["FROM"] = df.index.get_level_values("TECHNOLOGY").str[3:6]
        df["TO"] = df.index.get_level_values("TECHNOLOGY").str[8:11]
        df = df.drop_duplicates(subset=["FROM", "TO"], keep=False)  # intercountry
    else:
        df["FROM"] = df.index.get_level_values("TECHNOLOGY").str[3:8]
        df["TO"] = df.index.get_level_values("TECHNOLOGY").str[8:13]

    if df.empty:
        return pd.DataFrame(columns=["FROM", "TO", "YEAR", "VALUE"]).set_index(
            ["FROM", "TO", "YEAR"]
        )

    return (
        df.reset_index()[["FROM", "TO", "YEAR", "VALUE"]]
        .groupby(["FROM", "TO", "YEAR"])
        .sum()
    )


def calc_pwr_capacity(
    total_capacity_annual: pd.DataFrame, country: bool
) -> pd.DataFrame:

    df = total_capacity_annual.copy()

    df = df[
        (df.index.get_level_values("TECHNOLOGY").str.startswith("PWR"))
        & ~(df.index.get_level_values("TECHNOLOGY").str.contains("TRN"))
    ]

    df["TECH"] = df.index.get_level_values("TECHNOLOGY").str[3:6]

    if country:
        r = "COUNTRY"
        df[r] = df.index.get_level_values("TECHNOLOGY").str[6:9]
    else:
        r = "NODE"
        df[r] = df.index.get_level_values("TECHNOLOGY").str[6:11]

    return (
        df.reset_index()[["TECH", r, "YEAR", "VALUE"]]
        .groupby(["TECH", r, "YEAR"])
        .sum()
        .rename(columns={"TECH": "TECHNOLOGY"})
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        total_capacity_csv = snakemake.input.total_capacity
        pwr_node_save = snakemake.output.power_capacity_node
        trn_node_save = snakemake.output.transmission_capacity_node
        pwr_country_save = snakemake.output.power_capacity_country
        trn_country_save = snakemake.output.transmission_capacity_country
    else:
        total_capacity_csv = "results/India/results/TotalCapacityAnnual.csv"
        pwr_node_save = ""
        trn_node_save = ""
        pwr_country_save = ""
        trn_country_save = ""

    total_capacity_annual = pd.read_csv(total_capacity_csv, index_col=[0, 1, 2])

    pwr_node_capacity = calc_pwr_capacity(total_capacity_annual, country=False)
    trn_node_capacity = calc_trn_capacity(total_capacity_annual, country=False)
    pwr_country_capacity = calc_pwr_capacity(total_capacity_annual, country=True)
    trn_country_capacity = calc_trn_capacity(total_capacity_annual, country=True)

    pwr_node_capacity.to_csv(pwr_node_save, index=True)
    trn_node_capacity.to_csv(trn_node_save, index=True)
    pwr_country_capacity.to_csv(pwr_country_save, index=True)
    trn_country_capacity.to_csv(trn_country_save, index=True)
