"""Creates demand projection figures"""

import pandas as pd
import matplotlib.pyplot as plt

from regression import perform_regression
from data import get_historical_urban_pop_wb
from read import import_ember_elec, import_plexos_2015


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


if __name__ == "__main__":

    if "snakemake" in globals():
        file_plexos = snamkemake.inputs.plexos
        file_ember = snamkemake.inputs.ember
        regression_plot = snakemake.outputs.regression
    else:
        file_plexos = "resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx"
        file_ember = "resources/data/ember_yearly_electricity_data.csv"
        regression_plot = "regression.png"

    plexos = import_plexos_2015(file_plexos)
    ember = import_ember_elec(file_ember)
    urban = get_historical_urban_pop_wb(long=True)

    df = perform_regression(plexos, ember, urban)

    fig, axs = create_regression_plot(df, True)

    fig.savefig(regression_plot)
