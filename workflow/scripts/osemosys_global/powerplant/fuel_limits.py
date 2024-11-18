"""Applies fuel limits to mining technologies"""

import pandas as pd
from typing import Optional


"""

def apply_fuel_limits(region, years, output_data_dir, input_dir, max_fuel):

    mf_df = pd.read_csv(os.path.join(input_dir, "data", "fuel_limits.csv"))
    mf_df["TECHNOLOGY"] = "MIN" + mf_df["FUEL"] + mf_df["COUNTRY"]
    mf_df = mf_df[["TECHNOLOGY", "YEAR", "VALUE"]]

    tech_list = mf_df["TECHNOLOGY"].unique()
    mf_df_final = pd.DataFrame(
        list(itertools.product(tech_list, years)), columns=["TECHNOLOGY", "YEAR"]
    )
    mf_df_final = pd.merge(mf_df_final, mf_df, how="left", on=["TECHNOLOGY", "YEAR"])
    mf_df_final["VALUE"] = mf_df_final["VALUE"].astype(float)
    for each_tech in tech_list:
        mf_df_final.loc[mf_df_final["TECHNOLOGY"].isin([each_tech]), "VALUE"] = (
            mf_df_final.loc[mf_df_final["TECHNOLOGY"].isin([each_tech]), "VALUE"]
            .interpolate()
            .round(0)
        )

    mf_df_final["REGION"] = region
    mf_df_final = mf_df_final[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    mf_df_final.dropna(inplace=True)
    mf_df_final.to_csv(
        os.path.join(output_data_dir, "TotalTechnologyAnnualActivityUpperLimit.csv"),
        index=None,
    )
"""


def import_fuel_limits(f: str) -> pd.DataFrame:
    """Imports contry level fuel limits.

    fuel_limits.csv
    """
    return pd.read_csv(f)


def import_set(f: str) -> pd.Series:
    s = pd.read_csv(f).squeeze()
    if isinstance(s, pd.Series):
        return s
    else:
        # single element series
        return pd.Series(s)


def get_template_fuel_limit(
    regions: pd.Series, techs: pd.Series, years: pd.Series
) -> pd.Series:

    r = regions.unique().tolist()
    t = (
        techs[(techs.str.startswith("MIN")) & ~(techs.str.endswith("INT"))]
        .unique()
        .tolist()
    )
    y = years.unique().tolist()

    idx = pd.MultiIndex.from_product([r, t, y], names=["REGION", "TECHNOLOGY", "YEAR"])
    s = pd.Series(index=idx).fillna(0)
    s.name = "VALUE"

    return s


def get_user_fuel_limits(
    regions: pd.Series, fuel_limits: Optional[pd.DataFrame] = None
) -> pd.Series:

    r = regions.unique().tolist()
    assert len(r) == 1
    r = r[0]

    if fuel_limits.empty:
        s = pd.Series(index=["REGION", "TECHNOLOGY", "YEAR"])
        s.name = "VALUE"
        return s

    # format limits dataframe
    limits = fuel_limits.copy()
    limits = limits[limits.VALUE > 0.001].copy()  # assume these are zero
    limits["TECHNOLOGY"] = "MIN" + limits.FUEL + limits.COUNTRY
    limits = limits[["TECHNOLOGY", "YEAR", "VALUE"]].set_index(["TECHNOLOGY", "YEAR"])

    techs = limits.index.get_level_values("TECHNOLOGY").unique()

    dfs = []

    # need to interpolate each tehnology seperatly, as same years are not guaranteed
    for tech in techs:
        df = limits.loc[tech]
        idx = pd.Index(range(min(df.index), max(df.index) + 1))
        df = df.reindex(idx)
        df["VALUE"] = df.VALUE.interpolate(
            method="linear", limit_direction="both"
        ).round(1)
        df = df.reset_index().rename(columns={"index": "YEAR"})
        df["TECHNOLOGY"] = tech
        dfs.append(df[["TECHNOLOGY", "YEAR", "VALUE"]])

    final = pd.concat(dfs)

    final["REGION"] = r

    return final.set_index(["REGION", "TECHNOLOGY", "YEAR"]).squeeze()


def merge_template_user_limits(
    template: pd.Series, user_limits: pd.Series
) -> pd.Series:
    return user_limits.combine_first(template)


if __name__ == "__main__":

    if "snakemake" in globals():
        technology_csv = snakemake.input.technology_csv
        fuel_limit_csv = snakemake.input.fuel_limit_csv
        region_csv = snakemake.input.region_csv
        years_csv = snakemake.input.year_csv
        activity_upper_limit_csv = snakemake.output.annual_limit_csv
    else:
        technology_csv = "results/India/data/TECHNOLOGY.csv"
        fuel_limit_csv = "resources/data/fuel_limits.csv"
        region_csv = "results/India/data/REGION.csv"
        years_csv = "results/India/data/YEAR.csv"
        activity_upper_limit_csv = ""

    regions = import_set(region_csv)
    techs = import_set(technology_csv)
    years = import_set(years_csv)
    fuel_limits = import_fuel_limits(fuel_limit_csv)

    template = get_template_fuel_limit(regions, techs, years)
    user_limits = get_user_fuel_limits(regions, fuel_limits)

    activity_upper_limit = merge_template_user_limits(template, user_limits)

    activity_upper_limit.to_csv(activity_upper_limit_csv, index=True)
