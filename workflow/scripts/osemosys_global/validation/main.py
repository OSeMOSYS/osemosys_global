"""Main validation module 

For each county in the model scope, a figure is generated that will plot 
capacity, generation, and emissions for any year in 2022 and earlier. 
"""

import pandas as pd
from pathlib import Path
from utils import plot_gen_cap, format_og_data
from functools import partial
import eia
import irena
import ember

import logging

logger = logging.getLogger(__name__)

###
# getters
###


def get_generation_funcs(datasource: str) -> dict[str, callable]:
    match datasource:
        case "eia" | "EIA" | "Eia":
            return {
                "getter": eia.get_eia_generation,
                "formatter": partial(format_og_data, mapper=eia.OG_GEN_NAME_MAPPER),
                "plotter": plot_gen_cap,
            }
        case "irena" | "IRENA" | "Irena":
            return {
                "getter": irena.get_irena_generation,
                "formatter": partial(format_og_data, mapper=irena.OG_NAME_MAPPER),
                "plotter": plot_gen_cap,
            }
        case "ember" | "EMBER" | "Ember":
            return {
                "getter": ember.get_ember_generation,
                "formatter": partial(format_og_data, mapper=ember.OG_NAME_MAPPER),
                "plotter": plot_gen_cap,
            }
        case _:
            raise KeyError


def get_capacity_funcs(datasource: str) -> dict[str, callable]:
    match datasource:
        case "eia" | "EIA" | "Eia":
            return {
                "getter": eia.get_eia_capacity,
                "formatter": partial(format_og_data, mapper=eia.OG_CAP_NAME_MAPPER),
                "plotter": plot_gen_cap,
            }
        case "irena" | "IRENA" | "Irena":
            return {
                "getter": irena.get_irena_capacity,
                "formatter": partial(format_og_data, mapper=irena.OG_NAME_MAPPER),
                "plotter": plot_gen_cap,
            }
        case "ember" | "EMBER" | "Ember":
            return {
                "getter": ember.get_ember_capacity,
                "formatter": partial(format_og_data, mapper=ember.OG_NAME_MAPPER),
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
        datasource = "eia"
        variable = "capacity"
        result_dir = "results/India/results"
        data_file = "resources/data/validation/eia_capacity.json"
        options = {}
        # options = {"iso_codes": "resources/data/validation/iso.csv"}

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
        gen = funcs["plotter"](modelled, actual, variable, datasource)
        for country, (fig, _) in gen.items():
            p = Path(validation_results, country, variable)
            if not p.exists():
                p.mkdir(parents=True)
            f = Path(p, f"{datasource}.png")
            fig.tight_layout()
            fig.savefig(str(f))
