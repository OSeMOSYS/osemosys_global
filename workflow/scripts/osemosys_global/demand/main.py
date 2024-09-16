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
        file_plexos = snamkemake.inputs.plexos
        file_plexos_demand = snakemake.inputs.plexos_demand
        file_iamc_gdp = snamkemake.inputs.iamc_gdp
        file_iamc_pop = snamkemake.inputs.iamc_pop
        file_iamc_urb = snamkemake.inputs.iamc_urb
        file_iamc_missing = snamkemake.inputs.iamc_missing
        file_td_losses = snamkemake.inputs.td_losses
        file_ember = snamkemake.inputs.ember
        save_csv = snakemake.outputs.save
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
        save_csv = "SpecifiedAnnualDemand.csv"

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

    df.to_csv(save_csv, index=False)
