"""Applies fuel limits to mining technologies"""

import pandas as pd
from typing import Optional


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


def get_template_fuel_limit(regions: pd.Series, techs: pd.Series) -> pd.Series:

    r = regions.unique().tolist()
    t = (
        techs[(techs.str.startswith("MIN")) & ~(techs.str.endswith("INT"))]
        .unique()
        .tolist()
    )

    idx = pd.MultiIndex.from_product([r, t], names=["REGION", "TECHNOLOGY"])
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
        s = pd.Series(index=["REGION", "TECHNOLOGY"])
        s.name = "VALUE"
        return s

    limits = fuel_limits.copy()
    limits = limits[limits.VALUE > 0.001].copy()  # assume these are zero

    # not sure what the year column represents?
    limits = limits.drop(columns="YEAR").drop_duplicates(
        subset=["FUEL", "COUNTRY"], keep="last"
    )

    limits["TECHNOLOGY"] = "MIN" + limits.FUEL + limits.COUNTRY

    s = limits[["TECHNOLOGY", "VALUE"]].copy()
    s["REGION"] = r

    return s.set_index(["REGION", "TECHNOLOGY"]).squeeze()


def merge_template_user_limits(
    template: pd.Series, user_limits: pd.Series
) -> pd.Series:
    return user_limits.combine_first(template)


if __name__ == "__main__":

    if "snakemake" in globals():
        technology_csv = snakemake.input.technology_csv
        fuel_limit_csv = snakemake.input.fuel_limit_csv
        region_csv = snakemake.input.region_csv
        model_period_limit_csv = snakemake.output.model_period_limit_csv
    else:
        technology_csv = "results/India/data/TECHNOLOGY.csv"
        fuel_limit_csv = "resources/data/fuel_limits.csv"
        region_csv = "results/India/data/REGION.csv"
        model_period_limit_csv = ""

    regions = import_set(region_csv)
    techs = import_set(technology_csv)
    fuel_limits = import_fuel_limits(fuel_limit_csv)

    template = get_template_fuel_limit(regions, techs)
    user_limits = get_user_fuel_limits(regions, fuel_limits)

    model_period_limit = merge_template_user_limits(template, user_limits)

    model_period_limit.to_csv(model_period_limit_csv, index=True)
