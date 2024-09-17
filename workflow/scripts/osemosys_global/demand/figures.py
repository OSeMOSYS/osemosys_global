"""Creates demand projection figures"""

import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional
from datetime import datetime

from regression import perform_regression, get_regression_coefficients
from projection import perform_country_projection_step, _get_base_data
from data import get_historical_urban_pop_wb, get_iamc_data
from read import import_ember_elec, import_plexos_2015, import_iamc, import_iamc_missing
from constants import SPATIAL_RESOLUTION
from spatial import get_spatial_mapping_country


def create_regression_plot(df: pd.DataFrame, urbanization: bool) -> tuple:
    """Creates regression plot for continents

    df: pd.DataFrame
        regressed dataframe
    """

    num_regions = len(df.index.unique())

    fig, axs = plt.subplots(num_regions, figsize=(8, 6 * num_regions))

    for i, region in enumerate(df.index.unique()):

        data = df.loc[region]

        b = data["WB_GDPppp"]
        c = data["ember_Elec"]
        x = data.loc[region, "WB_GDPppp"]
        m = data.loc[region, "coef_GDPppp"]
        z = data.loc[region, "intercept"]

        data.plot.scatter(
            "WB_GDPppp",
            "ember_Elec",
            color="blue",
            alpha=0.5,
            label="GDPppp per capita",
            ax=axs[i],
        )

        if urbanization:

            ax = axs[i].twiny()
            w = data.loc[region, "R2_GDPppp_Urb/Elec"].unique()
            data.plot.scatter(
                "WB_Urb",
                "ember_Elec",
                color="red",
                alpha=0.5,
                label="Urban population",
                ax=ax,
            )
            ax.set_xlabel("Urban population (% of total population)")
            axs[i].legend(bbox_to_anchor=(1, 1))

        else:

            w = data.loc[region, "R2_GDPppp/Elec"].unique()
            data.plot(x, m * x + z, color="black", alpha=0.5, ax=axs[i])
            plt.title(f"{region} R2 = {w}")

        axs[i].legend(bbox_to_anchor=(0.31, 0.9))
        axs[i].set_title(f"{region} R2 = {w}")
        axs[i].set_ylabel("Electricity demand per capita (kWh)")
        axs[i].set_xlabel("GDPppp per capita (constant 2017 international $)")

    fig.tight_layout()

    return fig, axs


def create_demand_plot(
    plexos: pd.DataFrame,
    base: pd.DataFrame,
    reg: pd.DataFrame,
    dem: pd.DataFrame,
    urbanization: Optional[bool] = True,
) -> tuple:
    """Creates demand plot for continents

    df: pd.DataFrame
        demand dataframe
    """

    spatial_mapping = get_spatial_mapping_country(plexos)

    num_regions = len(spatial_mapping[SPATIAL_RESOLUTION].unique())

    fig, axs = plt.subplots(num_regions, figsize=(8, 5 * num_regions))

    for i, region in enumerate(spatial_mapping[SPATIAL_RESOLUTION].unique()):

        ctry_base = base[base[SPATIAL_RESOLUTION] == region]
        ctry_dem = dem[dem[SPATIAL_RESOLUTION] == region]
        ctry_reg = reg.loc[region]

        ctry_reg.plot.scatter(
            "WB_GDPppp",
            "ember_Elec",
            color="grey",
            alpha=0.5,
            label=f"2000 - {datetime.now().year}",
            ax=axs[i],
        )

        if not urbanization:
            pass
            # x = ctry_reg["WB_GDPppp"]
            # m = ctry_reg["coef_GDPppp"].unique()
            # b = ctry_reg["intercept"].unique()
            # plot(x, m * x + b, color="black", alpha=0.5)

        years = [2035, 2050, 2100]
        colours = ["green", "blue", "red"]

        for year, colour in zip(years, colours):

            temp = pd.concat(
                [
                    ctry_base[year].to_frame(name="base"),
                    ctry_dem[year].to_frame(name="dem"),
                ],
                axis=1,
            )

            temp.plot.scatter(
                "base",
                "dem",
                color=colour,
                alpha=0.5,
                label=year,
                ax=axs[i],
            )

            if not urbanization:
                # x_new = ctry_base[year]
                # plot(x_new, m * x + b, color="black", alpha=0.5)
                pass

        axs[i].set_title(f"Demand projection {region}")
        axs[i].set_ylabel("Electricity demand per capita (kWh)")
        axs[i].set_xlabel("GDPppp per capita")

    fig.tight_layout()

    return fig, axs


if __name__ == "__main__":

    if "snakemake" in globals():
        file_plexos = snakemake.input.plexos
        file_ember = snakemake.input.ember
        file_iamc_gdp = snakemake.input.iamc_gdp
        file_iamc_pop = snakemake.input.iamc_pop
        file_iamc_urb = snakemake.input.iamc_urb
        file_iamc_missing = snakemake.input.iamc_missing
        regression_plot = snakemake.output.regression
        projection_plot = snakemake.output.projection
    else:
        file_plexos = "resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx"
        file_ember = "resources/data/ember_yearly_electricity_data.csv"
        file_iamc_gdp = "resources/data/iamc_db_GDPppp_Countries.xlsx"
        file_iamc_pop = "resources/data/iamc_db_POP_Countries.xlsx"
        file_iamc_urb = "resources/data/iamc_db_URB_Countries.xlsx"
        file_iamc_missing = (
            "resources/data/iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx"
        )
        regression_plot = "regression.png"
        projection_plot = "projection.png"

    plexos = import_plexos_2015(file_plexos)
    ember = import_ember_elec(file_ember)
    urban = get_historical_urban_pop_wb(long=True)

    iamc_gdp_orig = import_iamc(file_iamc_gdp)
    iamc_pop_orig = import_iamc(file_iamc_pop)
    iamc_urb_orig = import_iamc(file_iamc_urb)
    iamc_gdp_missing = import_iamc_missing(file_iamc_missing, "gdp")
    iamc_pop_missing = import_iamc_missing(file_iamc_missing, "pop")
    iamc_urb_missing = import_iamc_missing(file_iamc_missing, "urb")

    iamc_gdp = get_iamc_data(plexos, iamc_gdp_orig, iamc_gdp_missing, "gdp")
    iamc_pop = get_iamc_data(plexos, iamc_pop_orig, iamc_pop_missing, "pop")
    iamc_urb = get_iamc_data(plexos, iamc_urb_orig, iamc_urb_missing, "urb")

    # create regression plots

    reg = perform_regression(plexos, ember, urban)

    fig, axs = create_regression_plot(reg, True)

    fig.savefig(regression_plot)

    # create demand projection plots

    dem = perform_country_projection_step(
        reg, iamc_gdp, iamc_pop, iamc_urb, per_capita=True
    )
    lr_coef = get_regression_coefficients(reg, True)
    dem_base = _get_base_data(iamc_gdp, iamc_pop, lr_coef)

    fig, axs = create_demand_plot(plexos, dem_base, reg, dem)

    fig.savefig(projection_plot)
