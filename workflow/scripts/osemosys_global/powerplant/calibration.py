import pandas as pd

def apply_calibration(calibration, oar_df, accumulated_annual_demand_df, fuel_set):

    cal_df = oar_df.loc[
        (oar_df["TECHNOLOGY"].str.startswith("PWR"))
        & ~(oar_df["TECHNOLOGY"].str.startswith("PWRTRN"))
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

        cal_dem_df = cal_df[["REGION", "FUEL", "YEAR", "DEMAND"]]
        cal_dem_final = cal_dem_df.rename(columns={"DEMAND": "VALUE"})
        cal_dem_final.drop_duplicates(inplace=True)
        cal_dem_final.dropna(inplace=True)
        
        accumulated_annual_demand_df = pd.concat([accumulated_annual_demand_df, 
                                                  cal_dem_final])

        cal_fuels = list(cal_dem_final["FUEL"].unique())
        fuel_set = list(fuel_set["VALUE"].unique()) + cal_fuels
        fuel_set = pd.DataFrame(fuel_set, columns=["VALUE"])
        
        return fuel_set, oar_df, accumulated_annual_demand_df