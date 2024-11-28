"""Creates demand projections"""

import pandas as pd
from read import (
    import_ember_elec,
    import_hourly_demand,
    import_iamc,
    import_iamc_missing,
    import_plexos_2015,
    import_td_losses,
)
from data import get_historical_urban_pop_wb, get_iamc_data, format_for_writing
from regression import perform_regression
from projection import perform_node_projections
from custom import (
    get_custom_demand_data,
    import_custom_demand_data,
    merge_default_custom_data,
)


def main(
    plexos: pd.DataFrame,
    ember: pd.DataFrame,
    plexos_demand: pd.DataFrame,
    iamc_gdp: pd.DataFrame,
    iamc_pop: pd.DataFrame,
    iamc_urb: pd.DataFrame,
    td_losses: pd.DataFrame,
) -> pd.DataFrame:

    wb_urban = get_historical_urban_pop_wb(long=True)

    regression = perform_regression(plexos, ember, wb_urban)

    projection = perform_node_projections(
        lr=regression,
        plexos=plexos,
        plexos_demand=plexos_demand,
        iamc_gdp=iamc_gdp,
        iamc_pop=iamc_pop,
        iamc_urb=iamc_urb,
        td_losses=td_losses,
    )

    df = format_for_writing(projection)

    return df


if __name__ == "__main__":

    # gets file paths
    if "snakemake" in globals():
        file_plexos = snakemake.input.plexos
        file_plexos_demand = snakemake.input.plexos_demand
        file_iamc_gdp = snakemake.input.iamc_gdp
        file_iamc_pop = snakemake.input.iamc_pop
        file_iamc_urb = snakemake.input.iamc_urb
        file_iamc_missing = snakemake.input.iamc_missing
        file_td_losses = snakemake.input.td_losses
        file_ember = snakemake.input.ember
        csv = snakemake.output.csv_files
        start_year = snakemake.params.start_year
        end_year = snakemake.params.end_year
        custom_nodes = snakemake.params.custom_nodes
        custom_nodes_data = snakemake.input.custom_nodes
    else:
        file_plexos = "resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx"
        file_plexos_demand = "resources/data/All_Demand_UTC_2015.csv"
        file_iamc_gdp = "resources/data/iamc_db_GDPppp_Countries.xlsx"
        file_iamc_pop = "resources/data/iamc_db_POP_Countries.xlsx"
        file_iamc_urb = "resources/data/iamc_db_URB_Countries.xlsx"
        file_iamc_missing = (
            "resources/data/iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx"
        )
        file_td_losses = "resources/data/T&D Losses.xlsx"
        file_ember = "resources/data/ember_yearly_electricity_data.csv"
        csv = "SpecifiedAnnualDemand.csv"
        start_year = 2020
        end_year = 2050
        custom_nodes = ["INDWE", "INDEA", "INDNE", "INDNO", "INDSO"]
        custom_nodes_data = "resources/data/custom_nodes/specified_annual_demand.csv"

    # first bring together original and missing iamc data
    plexos = import_plexos_2015(file_plexos)

    iamc_gdp_orig = import_iamc(file_iamc_gdp)
    iamc_pop_orig = import_iamc(file_iamc_pop)
    iamc_urb_orig = import_iamc(file_iamc_urb)
    iamc_gdp_missing = import_iamc_missing(file_iamc_missing, "gdp")
    iamc_pop_missing = import_iamc_missing(file_iamc_missing, "pop")
    iamc_urb_missing = import_iamc_missing(file_iamc_missing, "urb")

    iamc_gdp = get_iamc_data(plexos, iamc_gdp_orig, iamc_gdp_missing, "gdp")
    iamc_pop = get_iamc_data(plexos, iamc_pop_orig, iamc_pop_missing, "pop")
    iamc_urb = get_iamc_data(plexos, iamc_urb_orig, iamc_urb_missing, "urb")

    # perfrom regression and projection
    input_data = {
        "plexos": plexos,
        "ember": import_ember_elec(file_ember),
        "plexos_demand": import_hourly_demand(file_plexos_demand),
        "iamc_gdp": iamc_gdp,
        "iamc_pop": iamc_pop,
        "iamc_urb": iamc_urb,
        "td_losses": import_td_losses(file_td_losses),
    }

    df = main(**input_data)

    df = df[df.YEAR.isin(range(start_year, end_year + 1))]

    all_custom = import_custom_demand_data(custom_nodes_data)
    custom = get_custom_demand_data(all_custom, start_year, end_year)
    df = merge_default_custom_data(df, custom)

    df.to_csv(csv, index=False)
