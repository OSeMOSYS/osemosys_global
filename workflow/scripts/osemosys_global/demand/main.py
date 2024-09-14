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
from regression import perform_regression


def main(
    plexos: pd.DataFrame,
    ember: pd.DataFrame,
    demand: pd.DataFrame,
    iamc_gdp: pd.DataFrame,
    iamc_pop: pd.DataFrame,
    iamc_urb: pd.DataFrame,
    iamc_missing: pd.DataFrame,
    td_losses: pd.DataFrame,
) -> pd.DataFrame:
    
    regression = perform_regression(plexos, ember)
    
    gdp = 


if __name__ == "__main__":

    if "snakemake" in globals():
        file_plexos = snamkemake.inputs.plexos
        file_plexos_demand = snakemake.inputs.plexos_demand
        file_iamc_gdp = snamkemake.inputs.iamc_gdp
        file_iamc_pop = snamkemake.inputs.iamc_pop
        file_iamc_urb = snamkemake.inputs.iamc_urb
        file_iamc_missing = snamkemake.inputs.iamc_missing
        file_td_losses = snamkemake.inputs.td_losses
        file_ember = snamkemake.inputs.ember
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

    input_data = {
        "plexos": import_plexos_2015(file_plexos),
        "ember": import_ember_elec(file_ember),
        "demand": import_hourly_demand(file_plexos_demand),
        "iamc_gdp": import_iamc(file_iamc_gdp),
        "iamc_pop": import_iamc(file_iamc_pop),
        "iamc_urb": import_iamc(file_iamc_urb),
        "iamc_missing": import_iamc_missing(file_iamc_missing),
        "td_losses": import_td_losses(file_td_losses),
    }

    df = main(**input_data)
