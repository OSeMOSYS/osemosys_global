"""Main validation module 

For each county in the model scope, a figure is generated that will plot 
capacity, generation, and emissions for any year in 2022 and earlier. 
"""

import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional
from pathlib import Path
import eia
import irena
import ember

import logging

logger = logging.getLogger(__name__)

###
# plotters
###


def plot_gen_cap(
    modelled: pd.DataFrame,
    actual: pd.DataFrame,
    variable: str,
    dataset_name: Optional[str] = None,
) -> dict[str, tuple[plt.figure, plt.axes]]:

    def _join_data(
        modelled: pd.DataFrame, actual: pd.DataFrame, dataset_name: Optional[str] = None
    ) -> pd.DataFrame:

        if not dataset_name:
            dataset_name = "ACTUAL"

        modelled = modelled.rename(columns={"VALUE": "OSeMOSYS"})
        actual = actual.rename(columns={"VALUE": dataset_name})
        df = modelled.join(actual)

        assert len(df.index.get_level_values("REGION").unique()) == 1

        return df.droplevel("REGION")

    assert modelled.index.names == actual.index.names

    if variable == "generation":
        units = "PJ"
    elif variable == "capacity":
        units = "GW"
    else:
        raise ValueError(
            f"Variable must be one of ['generation', 'capacity']. Recieved {variable}"
        )

    df = _join_data(modelled, actual, dataset_name).reset_index()
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
            title = f"{country} {variable.capitalize()} in {year}"
            if n_rows > 1:
                ax = axs[i]
            else:
                ax = axs
            df_year.plot(
                kind="bar", ax=ax, rot=45, title=title, xlabel="", ylabel=units
            )

        data[country] = (fig, axs)

    return data


###
# getters
###


def get_generation_funcs(datasource: str) -> dict[str, callable]:
    match datasource:
        case "eia" | "EIA" | "Eia":
            return {
                "getter": eia.get_eia_generation,
                "formatter": eia.format_og_generation,
                "plotter": plot_gen_cap,
            }
        case "irena" | "IRENA" | "Irena":
            return {
                "getter": irena.get_irena_generation,
                "formatter": irena.format_og_data,
                "plotter": plot_gen_cap,
            }
        case "ember" | "EMBER" | "Ember":
            return {
                "getter": ember.get_ember_generation,
                "formatter": ember.format_og_data,
                "plotter": plot_gen_cap,
            }
        case _:
            raise KeyError


def get_capacity_funcs(datasource: str) -> dict[str, callable]:
    match datasource:
        case "eia" | "EIA" | "Eia":
            return {
                "getter": eia.get_eia_capacity,
                "formatter": eia.format_og_capacity,
                "plotter": plot_gen_cap,
            }
        case "irena" | "IRENA" | "Irena":
            return {
                "getter": irena.get_irena_capacity,
                "formatter": irena.format_og_data,
                "plotter": plot_gen_cap,
            }
        case "ember" | "EMBER" | "Ember":
            return {
                "getter": ember.get_ember_capacity,
                "formatter": ember.format_og_data,
                "plotter": plot_gen_cap,
            }
        case _:
            raise KeyError


###
# entry point
###

if __name__ == "__main__":
    if "snakemake" in globals():
        raise NotImplementedError
    else:
        datasource = "irena"
        variable = "capacity"
        result_dir = "results/India/results"
        data_file = "resources/data/validation/irena_capacity.csv"
        options = {}
        options = {"iso_codes": "resources/data/validation/iso.csv"}

    csv_results = Path(result_dir)
    validation_results = Path(csv_results, "..", "validation")

    # get data based on validation variable and datasource

    if variable == "generation":
        og_result = "ProductionByTechnologyAnnual"
        funcs = get_generation_funcs(datasource)
    elif variable == "capacity":
        og_result = "TotalCapacityAnnual"
        funcs = get_capacity_funcs(datasource)
    else:
        raise NotImplementedError

    # perform the validation

    try:
        actual = funcs["getter"](data_file, **options)
        modelled = pd.read_csv(Path(result_dir, f"{og_result}.csv"))
        modelled = funcs["formatter"](modelled)
    except KeyError as e:
        actual = None
        modelled = None
        logger.error(f"No validation for {variable} from {datasource}: \n{e}")

    if isinstance(actual, pd.DataFrame) and isinstance(modelled, pd.DataFrame):
        gen = funcs["plotter"](modelled, actual, datasource)
        for country, (fig, _) in gen.items():
            p = Path(validation_results, country, variable)
            if not p.exists():
                p.mkdir(parents=True)
            f = Path(p, f"{datasource}.png")
            fig.tight_layout()
            fig.savefig(str(f))
