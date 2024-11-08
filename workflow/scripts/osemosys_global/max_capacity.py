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

    input_dir = config_paths.input_dir
    output_data_dir = config_paths.output_data_dir
    region = config.region_name
    years = config.get_years()
    remove_nodes = config.get("nodes_to_remove")
    max_build = config.get("powerplant_build_rates")
    max_fuel = config.get("fuel_limits")
    calibration = config.get("calibration")
    re_targets = config.get("re_targets")
    
    apply_build_rates(region, years, output_data_dir, input_dir, max_build)
    apply_fuel_limits(region, years, output_data_dir, input_dir, max_fuel)
    apply_calibration(region, years, output_data_dir, calibration)
    apply_re_targets(region, years, output_data_dir, re_targets, remove_nodes)

def apply_build_rates(region, years, output_data_dir, input_dir, max_build):
    max_build_df = pd.read_csv(
        os.path.join(input_dir, "data", "powerplant_build_rates.csv")
    )
    max_build_df = max_build_df.loc[
        max_build_df.index.repeat(
            max_build_df["END_YEAR"] + 1 - max_build_df["START_YEAR"]
        )
    ]
    max_build_df["YEAR"] = (
        max_build_df.groupby(level=0).cumcount() + max_build_df["START_YEAR"]
    )
    max_build_df = max_build_df.reset_index(drop=True)
    max_build_df = max_build_df[["TYPE", "METHOD", "MAX_BUILD", "YEAR", "COUNTRY"]]
    max_build_df["TYPE"] = max_build_df["TYPE"].str[0:3]
    # Create a list of powerplant technologies
    tech_set = pd.read_csv(os.path.join(output_data_dir, "TECHNOLOGY.csv"))
    pwr_tech_list = [x for x in list(tech_set["VALUE"]) if x.startswith("PWR")]

    # Create scaffold dataframe with all powerplant technologies for all years
    df_techs = pd.DataFrame(
        list(itertools.product(pwr_tech_list, years)), columns=["TECHNOLOGY", "YEAR"]
    )

    # Filter out technologies for which a max. capacity investment has already
    # been set
    df_max_cap_inv = pd.read_csv(
        os.path.join(output_data_dir, "TotalAnnualMaxCapacityInvestment.csv")
    )
    max_cap_inv_techs = list(df_max_cap_inv["TECHNOLOGY"].unique())
    df_techs = df_techs[~(df_techs["TECHNOLOGY"].isin(max_cap_inv_techs))]

    # Create dataframe of max capacity by technology
    if os.path.isfile(os.path.join(output_data_dir, "TotalAnnualMaxCapacity.csv")):
        df_max_cap = pd.read_csv(
            os.path.join(output_data_dir, "TotalAnnualMaxCapacity.csv")
        )
        df_max_cap = df_max_cap.loc[df_max_cap["TECHNOLOGY"].str.startswith("PWR")]
    else:
        df_max_cap = pd.DataFrame(columns=["REGION", "TECHNOLOGY", "YEAR", "VALUE"])
    df_max_cap = pd.merge(
        left=df_techs, right=df_max_cap, on=["TECHNOLOGY", "YEAR"], how="left"
    )
    df_max_cap["REGION"] = region
    df_max_cap["TYPE"] = df_max_cap["TECHNOLOGY"].str[3:6]
    df_max_cap["COUNTRY"] = df_max_cap["TECHNOLOGY"].str[6:9]
    df_max_cap = pd.merge(
        left=df_max_cap, right=max_build_df, on=["TYPE", "COUNTRY", "YEAR"], how="left"
    )
    df_max_cap.loc[df_max_cap["METHOD"].isin(["ABS"]), "VALUE"] = df_max_cap[
        "MAX_BUILD"
    ]
    df_max_cap.dropna(inplace=True)
    df_max_cap.loc[df_max_cap["METHOD"].isin(["PCT"]), "VALUE"] = (
        df_max_cap["VALUE"] * df_max_cap["MAX_BUILD"] / 100
    )

    df_max_cap = df_max_cap[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    df_max_cap = pd.concat([df_max_cap_inv, df_max_cap], ignore_index=True)
    df_max_cap["VALUE"] = df_max_cap["VALUE"].astype(float).round(3)
    df_max_cap.to_csv(
        os.path.join(output_data_dir, "TotalAnnualMaxCapacityInvestment.csv"),
        index=None,
    )


def apply_fuel_limits(region, years, output_data_dir, input_dir, max_fuel):

    mf_df = pd.read_csv(os.path.join(input_dir, "data", "fuel_limits.csv"))
    mf_df["TECHNOLOGY"] = "MIN" + mf_df["FUEL"] + mf_df["COUNTRY"]
    mf_df = mf_df[["TECHNOLOGY", "YEAR", "VALUE"]]

    tech_list = mf_df["TECHNOLOGY"].unique()
    mf_df_final = pd.DataFrame(
        list(itertools.product(tech_list, years)), columns=["TECHNOLOGY", "YEAR"]
    )
    mf_df_final = pd.merge(mf_df_final, mf_df, how="left", on=["TECHNOLOGY", "YEAR"])
    mf_df_final["VALUE"] = mf_df_final["VALUE"].astype(float)
    for each_tech in tech_list:
        mf_df_final.loc[mf_df_final["TECHNOLOGY"].isin([each_tech]), "VALUE"] = (
            mf_df_final.loc[mf_df_final["TECHNOLOGY"].isin([each_tech]), "VALUE"]
            .interpolate()
            .round(0)
        )

    mf_df_final["REGION"] = region
    mf_df_final = mf_df_final[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    mf_df_final.dropna(inplace=True)
    mf_df_final.to_csv(
        os.path.join(output_data_dir, "TotalTechnologyAnnualActivityUpperLimit.csv"),
        index=None,
    )

    # Model Period Activity Upper Limit for 'MINCOA***01'
    min_tech_df = pd.read_csv(os.path.join(output_data_dir, "TECHNOLOGY.csv"))
    min_tech = [
        x
        for x in min_tech_df["VALUE"].unique()
        if x.startswith("MINCOA")
        if x.endswith("01")
    ]
    min_tech_df_final = pd.DataFrame(columns=["REGION", "TECHNOLOGY", "VALUE"])
    min_tech_df_final["TECHNOLOGY"] = min_tech
    min_tech_df_final["REGION"] = region
    min_tech_df_final["VALUE"] = 0
    min_tech_df_final.to_csv(
        os.path.join(
            output_data_dir, "TotalTechnologyModelPeriodActivityUpperLimit.csv"
        ),
        index=None,
    )


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