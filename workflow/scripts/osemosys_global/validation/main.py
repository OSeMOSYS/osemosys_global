"""Main validation module 

For each county in the model scope, a figure is generated that will plot 
capacity, generation, and emissions for any year in 2022 and earlier. 

Capacity comparisons: 
- EIA International Energy Outlook (https://www.eia.gov/international/overview/world)

Generation comparisons: 
- EIA International Energy Outlook (https://www.eia.gov/international/overview/world)
"""

import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional
from pathlib import Path
import eia

import logging

logger = logging.getLogger(__name__)


def plot_generation(
    og: pd.DataFrame, real: pd.DataFrame, dataset_name: Optional[str] = None
) -> dict[str, tuple[plt.figure, plt.axes]]:

    def _join_data(
        og: pd.DataFrame, real: pd.DataFrame, dataset_name: Optional[str] = None
    ) -> pd.DataFrame:
        if not dataset_name:
            dataset_name = "ACTUAL"
        og = og.rename(columns={"VALUE": "OSeMOSYS"})
        real = real.rename(columns={"VALUE": dataset_name})
        df = og.join(real)
        assert len(df.index.get_level_values("REGION").unique()) == 1
        return df.droplevel("REGION")

    assert og.index.names == real.index.names

    df = _join_data(og, real, dataset_name).reset_index()
    df["TECH"] = df["TECHNOLOGY"].str[0:3]
    df["COUNTRY"] = df["TECHNOLOGY"].str[3:]

    data = {}

    countries = df.COUNTRY.unique()
    for country in countries:
        df_country = df[df.COUNTRY == country]
        years = df_country.YEAR.unique()
        n_rows = len(years)
        fig, axs = plt.subplots(n_rows, 1, figsize=(10, n_rows * 4))
        for i, year in enumerate(years):
            df_year = (
                df_country[df_country.YEAR == year]
                .drop(columns=["TECHNOLOGY", "YEAR", "COUNTRY"])
                .set_index("TECH")
            )
            title = f"{country} Generation in {year}"
            if n_rows > 1:
                ax = axs[i]
            else:
                ax = axs
            df_year.plot(kind="bar", ax=ax, rot=45, title=title, xlabel="", ylabel="PJ")

        data[country] = (fig, axs)

    return data


def get_generation_funcs(datasource: str) -> dict[str, callable]:
    match datasource:
        case "eia" | "EIA":
            return {
                "getter": eia.get_eia_generation,
                "formatter": eia.format_og_generation,
            }
        case _:
            raise KeyError


if __name__ == "__main__":
    if "snakemake" in globals():
        raise NotImplementedError
    else:
        datasource = "eia"
        result_dir = "results/India/results"
        capacity_file = "resources/data/validation/eia_capacity.json"
        generation_file = "resources/data/validation/eia_generation.json"

    csv_results = Path(result_dir)
    validation_results = Path(csv_results, "..", "validation")

    # generation validation

    try:
        funcs = get_generation_funcs(datasource)
        actual_gen = funcs["getter"](generation_file)
        og_gen = pd.read_csv(Path(result_dir, "ProductionByTechnologyAnnual.csv"))
        og_gen = funcs["formatter"](og_gen)
    except KeyError:
        actual_gen = None
        og_gen = None
        logger.info(f"No generation validation for {datasource}")

    if isinstance(actual_gen, pd.DataFrame) and isinstance(og_gen, pd.DataFrame):
        gen = plot_generation(og_gen, actual_gen, datasource)
        for country, (fig, _) in gen.items():
            p = Path(validation_results, country, "gen")
            if not p.exists():
                p.mkdir(parents=True)
            f = Path(p, f"{datasource}.png")
            fig.tight_layout()
            fig.savefig(str(f))
