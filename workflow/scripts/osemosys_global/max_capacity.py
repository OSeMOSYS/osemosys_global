import os
import pandas as pd
from configuration import ConfigFile, ConfigPaths
import itertools

# LOGGING
import logging
logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)


def main():
    """Creates capacity limits on renewable technologies."""

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    config = ConfigFile("config")

    output_data_dir = config_paths.output_data_dir
    region = config.region_name
    years = config.get_years()
    remove_nodes = config.get("nodes_to_remove")
    calibration = config.get("calibration")
    re_targets = config.get("re_targets")
    
    apply_calibration(region, years, output_data_dir, calibration)
    apply_re_targets(region, years, output_data_dir, re_targets, remove_nodes)

def apply_calibration(region, years, output_data_dir, calibration):

    oar_df = pd.read_csv(os.path.join(output_data_dir, "OutputActivityRatio.csv"))
    cal_df = oar_df.loc[
        (oar_df["TECHNOLOGY"].str.startswith("PWR"))
        & ~(oar_df["TECHNOLOGY"].str.startswith("PWRTRN"))
        & ~(oar_df["TECHNOLOGY"].str.startswith("PWRBAT"))
    ]

    if not calibration is None:
        for cal, cal_params in calibration.items():

            if cal[0:3] in ["GAS"]:
                cal_df.loc[
                    (
                        (cal_df["TECHNOLOGY"].str.startswith("PWROCG" + cal_params[1]))
                        | (
                            cal_df["TECHNOLOGY"].str.startswith(
                                "PWRCCG" + cal_params[1]
                            )
                        )
                    )
                    & (cal_df["YEAR"] == cal_params[2]),
                    "FUEL",
                ] = (
                    "CALGAS" + cal_params[1]
                )
                cal_df.loc[
                    (
                        (cal_df["TECHNOLOGY"].str.startswith("PWROCG" + cal_params[1]))
                        | (
                            cal_df["TECHNOLOGY"].str.startswith(
                                "PWRCCG" + cal_params[1]
                            )
                        )
                    )
                    & (cal_df["YEAR"] == cal_params[2]),
                    "DEMAND",
                ] = cal_params[0]
            else:
                cal_df.loc[
                    (
                        cal_df["TECHNOLOGY"].str.startswith(
                            "PWR" + cal[0:3] + cal_params[1]
                        )
                    )
                    & (cal_df["YEAR"] == cal_params[2]),
                    "FUEL",
                ] = (
                    "CAL" + cal[0:3] + cal_params[1]
                )
                cal_df.loc[
                    (
                        cal_df["TECHNOLOGY"].str.startswith(
                            "PWR" + cal[0:3] + cal_params[1]
                        )
                    )
                    & (cal_df["YEAR"] == cal_params[2]),
                    "DEMAND",
                ] = cal_params[0]

        cal_df = cal_df.loc[cal_df["FUEL"].str.startswith("CAL")]
        cal_df_final = cal_df[[x for x in cal_df.columns if x not in ["DEMAND"]]]
        oar_df = pd.concat([oar_df, cal_df_final])

        oar_df.to_csv(
            os.path.join(output_data_dir, "OutputActivityRatio.csv"), index=None
        )

        cal_dem_df = cal_df[["REGION", "FUEL", "YEAR", "DEMAND"]]
        cal_dem_final = cal_dem_df.rename(columns={"DEMAND": "VALUE"})
        cal_dem_final.drop_duplicates(inplace=True)
        cal_dem_final.dropna(inplace=True)
        cal_dem_final.to_csv(
            os.path.join(output_data_dir, "AccumulatedAnnualDemand.csv"), index=None
        )

        cal_fuels = list(cal_dem_final["FUEL"].unique())

        fuel_list_df = pd.read_csv(os.path.join(output_data_dir, "FUEL.csv"))
        fuel_list = list(fuel_list_df["VALUE"].unique()) + cal_fuels
        fuel_list_df_final = pd.DataFrame(fuel_list, columns=["VALUE"])
        fuel_list_df_final.to_csv(os.path.join(output_data_dir, "FUEL.csv"), index=None)
    else:
        cal_dem_final = pd.DataFrame(columns=["REGION", "FUEL", "YEAR", "VALUE"])
        cal_dem_final.to_csv(
            os.path.join(output_data_dir, "AccumulatedAnnualDemand.csv"), index=None
        )


def apply_re_targets(region, years, output_data_dir, re_targets, remove_nodes):
    """Apply Renewable Energy targets by country and year"""

    if not remove_nodes:
        remove_nodes = []

    oar_df = pd.read_csv(os.path.join(output_data_dir, "OutputActivityRatio.csv"))

    # List of RE technologies
    re_techs = ["BIO", "CSP", "GEO", "HYD", "SPV", "WON", "WOF", "WAS", "WAV"]
    re_df = oar_df.loc[
        (oar_df["TECHNOLOGY"].str.startswith("PWR"))
        & (oar_df["TECHNOLOGY"].str[3:6].isin(re_techs))
    ].copy()

    # Create dummy commodity starting with 'REN' for renewables
    re_df["FUEL"] = "REN" + re_df["FUEL"].str[3:6]
    oar_df = pd.concat([oar_df, re_df])
    oar_df.to_csv(os.path.join(output_data_dir, "OutputActivityRatio.csv"), index=None)

    # Create list of fuels
    re_fuels = list(re_df.loc[re_df["FUEL"].str.startswith("REN"), "FUEL"].unique())
    fuels_df = pd.read_csv(os.path.join(output_data_dir, "FUEL.csv"))
    fuels_ren_df = re_df[["FUEL"]].copy()
    fuels_ren_df.rename(columns={"FUEL": "VALUE"}, inplace=True)
    fuels_df = pd.concat([fuels_df, fuels_ren_df])
    fuels_df.drop_duplicates(inplace=True)
    fuels_df.to_csv(os.path.join(output_data_dir, "FUEL.csv"), index=None)

    # Create dataframe template to calculate SpecifiedAnnualDemand
    re_targets_df = pd.DataFrame(
        list(itertools.product(re_fuels, years)), columns=["FUEL", "YEAR"]
    )
    re_targets_df["COUNTRY"] = re_targets_df["FUEL"].str[3:6]
    re_targets_df = re_targets_df[["COUNTRY", "YEAR"]]

    if not re_targets is None:
        for re_params in re_targets:
            re_targets_df.loc[
                (re_targets_df["YEAR"].between(re_params[1], re_params[2]))
                & (re_targets_df["COUNTRY"].isin([re_params[0]])),
                "VALUE",
            ] = (
                re_params[3] / 100
            )

        re_targets_df = re_targets_df.pivot(
            index=["YEAR"], columns=["COUNTRY"], values="VALUE"
        ).reset_index()
        re_targets_df = re_targets_df.interpolate()

        # Drop all columns with only NaN
        re_targets_df.dropna(axis=1, how="all", inplace=True)

        # Melt to get 'COUNTRY' column and remove all rows with NaN
        re_targets_df = pd.melt(
            re_targets_df,
            id_vars=["YEAR"],
            value_vars=[x for x in re_targets_df.columns if x not in ["YEAR"]],
            var_name="COUNTRY",
            value_name="VALUE",
        )
        re_targets_df.dropna(axis=0, inplace=True)

        # Read 'SpecifiedAnnualDemand'
        sp_demand_df = pd.read_csv(
            os.path.join(output_data_dir, "SpecifiedAnnualDemand.csv")
        )
        sp_demand_df = sp_demand_df.loc[
            ~(sp_demand_df["FUEL"].str[3:8].isin(remove_nodes))
        ]
        sp_demand_df["COUNTRY"] = sp_demand_df["FUEL"].str[3:6]
        sp_demand_df = sp_demand_df.groupby(["YEAR", "COUNTRY"], as_index=False)[
            "VALUE"
        ].sum()
        sp_demand_df.rename(columns={"VALUE": "DEMAND"}, inplace=True)
        # Merge RE targets and SpecifiedAnnualDemand tables
        re_targets_df = pd.merge(
            re_targets_df, sp_demand_df, how="left", on=["COUNTRY", "YEAR"]
        )
        re_targets_df["VALUE"] = re_targets_df["VALUE"] * re_targets_df["DEMAND"]
        re_targets_df["FUEL"] = "REN" + re_targets_df["COUNTRY"]
        re_targets_df["REGION"] = region
        re_targets_df = re_targets_df[["REGION", "FUEL", "YEAR", "VALUE"]]

        # Read 'SpecifiedAnnualDemand'
        ac_demand_df = pd.read_csv(
            os.path.join(output_data_dir, "AccumulatedAnnualDemand.csv")
        )

        # Append RE target demands and write new AccumulatedAnnualDemand file
        ac_demand_df.drop_duplicates(inplace=True)
        ac_demand_df = pd.concat([ac_demand_df, re_targets_df])
        ac_demand_df.to_csv(
            os.path.join(output_data_dir, "AccumulatedAnnualDemand.csv"), index=None
        )


if __name__ == "__main__":
    main()
    logging.info("Max capacity limits sucessfully set")