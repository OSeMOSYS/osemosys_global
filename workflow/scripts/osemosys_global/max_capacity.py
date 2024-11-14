import os
import pandas as pd
from configuration import ConfigFile, ConfigPaths

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
    calibration = config.get("calibration")
    
    apply_calibration(region, years, output_data_dir, calibration)

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

if __name__ == "__main__":
    main()
    logging.info("Max capacity limits sucessfully set")