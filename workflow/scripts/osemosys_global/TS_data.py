import pandas as pd
import itertools
import seaborn as sns

sns.set()
import os

from osemosys_global.utils import apply_timeshift
from utils import apply_dtypes
from constants import SET_DTYPES
from datetime import datetime

def main(
        demand_df: pd.DataFrame,
        csp_df: pd.DataFrame,
        spv_df: pd.DataFrame,
        hyd_df: pd.DataFrame,
        won_df: pd.DataFrame,
        wof_df: pd.DataFrame,
        custom_sp_demand_profile: pd.DataFrame,
        csp_df_custom: pd.DataFrame,
        hyd_df_custom: pd.DataFrame,
        spv_df_custom: pd.DataFrame,
        wof_df_custom: pd.DataFrame,
        won_df_custom: pd.DataFrame,
        seasons_raw: dict,
        dayparts_raw: dict,
        ):
   
    seasonsData = []
    for s, months in seasons_raw.items():
        for month in months:
            seasonsData.append([month, s])
    seasons_df = pd.DataFrame(seasonsData, columns=["month", "season"])
    seasons_df = seasons_df.sort_values(by=["month"]).reset_index(drop=True)
    
    daypartData = []
    for dp, hr in dayparts_raw.items():
        daypartData.append([dp, hr[0], hr[1]])
    dayparts_df = pd.DataFrame(daypartData, columns=["daypart", "start_hour", "end_hour"])

    dayparts_df["start_hour"] = dayparts_df["start_hour"].map(
        lambda x: apply_timeshift(x, timeshift)
    )
    dayparts_df["end_hour"] = dayparts_df["end_hour"].map(
        lambda x: apply_timeshift(x, timeshift)
    )

    years = list(range(start_year, end_year + 1))
    
    # Read renewable profile files
    csp_df_custom.drop(["Datetime"], axis=1, inplace=True)
    csp_df = pd.concat([csp_df, csp_df_custom], axis=1)
    
    csp_df.name = "CSP"

    spv_df_custom.drop(["Datetime"], axis=1, inplace=True)
    spv_df = pd.concat([spv_df, spv_df_custom], axis=1)
    
    spv_df.name = "SPV"
    
    nodes = ["-".join(x.split("-")[1:]) for x in spv_df.columns if x not in ["Datetime"]]
    regions = [x for x in spv_df.columns if x not in ["Datetime"]]
    
    node_region_dict = dict(zip(nodes, regions))
        
    hyd_df = hyd_df.loc[hyd_df["NAME"].str.endswith("Capacity Scaler")]
    hyd_df.loc[:,"NAME"] = hyd_df["NAME"].str.split("_").str[0]

    hyd_df = pd.concat([hyd_df, hyd_df_custom])
    
    # Drop Brazil transmission nodes J1, J2, J3
    brazil_j_nodes = ["BRA-J1", "BRA-J2", "BRA-J3"]
    hyd_df = hyd_df.loc[~hyd_df["NAME"].isin(brazil_j_nodes)]
    hyd_df = hyd_df.set_index("NAME").T.reset_index()
    hyd_df.rename(columns={"index": "MONTH"}, inplace=True)
    
    hyd_df["MONTH"] = hyd_df["MONTH"].str.replace("M", "").astype(int)
    
    hyd_df_processed = pd.DataFrame(columns=["Datetime"])
    hyd_df_processed["Datetime"] = spv_df["Datetime"]
    hyd_df_processed["MONTH"] = (
        hyd_df_processed["Datetime"].str.split("/").str[1].astype(int)
    )
    hyd_df_processed = pd.merge(hyd_df_processed, hyd_df, how="left", on="MONTH")
    hyd_df_processed.drop(columns="MONTH", inplace=True)
    hyd_df_processed.rename(columns=node_region_dict, inplace=True)
    hyd_df_processed.name = "HYD"
    
    won_df_custom.drop(["Datetime"], axis=1, inplace=True)
    won_df = pd.concat([won_df, won_df_custom], axis=1)
    won_df.name = "WON"
    
    wof_df_custom.drop(["Datetime"], axis=1, inplace=True)
    wof_df = pd.concat([wof_df, wof_df_custom], axis=1)
    wof_df.name = "WOF"
       
    # ### Create 'output' directory if it doesn't exist
    if not os.path.exists(output_data_dir):
        os.makedirs(output_data_dir)
    
    
    def correct_datetime_formatting(time_str):
        """start hours are not padded with 00:00"""
        try:
            date_obj = datetime.strptime(time_str, "%d/%m/%Y %H:%M")
        except ValueError:
            try:
                date_obj = datetime.strptime(time_str, "%d/%m/%Y")
                date_obj = date_obj.replace(hour=0, minute=0)
            except ValueError as ex:
                print(ex)
                raise ValueError
        return date_obj
    
    
    # ### Create columns for year, month, day, hour, and day type
    
    # Convert datetime to year, month, day, and hour
    demand_df["Datetime"] = demand_df.Datetime.map(correct_datetime_formatting)
    demand_df["Datetime"] = pd.to_datetime(demand_df["Datetime"], format="%d/%m/%Y %H:%M")
    demand_df["Year"] = demand_df["Datetime"].dt.strftime("%Y").astype(int)
    demand_df["Month"] = demand_df["Datetime"].dt.strftime("%m").astype(int)
    demand_df["Day"] = demand_df["Datetime"].dt.strftime("%d").astype(int)
    demand_df["Hour"] = demand_df["Datetime"].dt.strftime("%H").astype(int)
    
    demand_nodes = [x for x in demand_df.columns if x != "Datetime"] + [
        y for y in custom_sp_demand_profile.iloc[:,3:].columns]
    
    demand_df = pd.merge(
        demand_df, custom_sp_demand_profile, how="left", on=["Month", "Day", "Hour"]
    )
    
    # Create column for weekday/weekend
    demand_df["Day-of-week"] = demand_df["Datetime"].dt.dayofweek.astype(str)
    demand_df.loc[demand_df["Day-of-week"].isin([0,1,2,3,4]), "Day-of-week"] = "WD"
    demand_df.loc[demand_df["Day-of-week"] != "WD", "Day-of-week"] = "WE"
    
    # ### Create dictionaries for 'seasons' and 'dayparts'
    
    seasons_dict = dict(zip(list(seasons_df["month"]), list(seasons_df["season"])))
    
    dayparts_dict = {
        i: [j, k]
        for i, j, k in zip(
            list(dayparts_df["daypart"]),
            list(dayparts_df["start_hour"]),
            list(dayparts_df["end_hour"]),
        )
    }
    
    
    # ### Create columns with 'seasons' and 'dayparts'
    
    
    demand_df["Season"] = demand_df["Month"]
    demand_df["Season"] = demand_df["Season"].replace(seasons_dict)
    
    demand_df["Hour"] = demand_df["Hour"].map(lambda x: apply_timeshift(int(x), timeshift))
    for daypart in dayparts_dict:
        if dayparts_dict[daypart][0] > dayparts_dict[daypart][1]:  # loops over 24hrs
            demand_df.loc[
                (demand_df["Hour"] >= dayparts_dict[daypart][0])
                | (demand_df["Hour"] < dayparts_dict[daypart][1]),
                "Daypart",
            ] = daypart
        else:
            demand_df.loc[
                (demand_df["Hour"] >= dayparts_dict[daypart][0])
                & (demand_df["Hour"] < dayparts_dict[daypart][1]),
                "Daypart",
            ] = daypart
    
    
    # ### Create column for timeslice with and without day-type
    
    
    if daytype:
        demand_df["TIMESLICE"] = (
            demand_df["Season"] + demand_df["Day-of-week"] + demand_df["Daypart"]
        )
    else:
        demand_df["TIMESLICE"] = demand_df["Season"] + demand_df["Daypart"]
    
    # ### Calculate YearSplit
    
    
    yearsplit = (
        demand_df["TIMESLICE"]
        .value_counts(normalize=True)
        .to_frame("VALUE")
        .round(4)
        .reset_index()
        .rename({"index": "TIMESLICE"}, axis=1)
    )
    
    yearsplit_final = pd.DataFrame(
        list(itertools.product(yearsplit["TIMESLICE"].unique(), years)),
        columns=["TIMESLICE", "YEAR"],
    )
    yearsplit_final = yearsplit_final.join(yearsplit.set_index("TIMESLICE"), on="TIMESLICE")
    yearsplit_final = apply_dtypes(yearsplit_final, "Year Split")
    yearsplit_final.to_csv(os.path.join(output_data_dir, "YearSplit.csv"), index=None)
    
    
    #  Calculate SpecifiedAnnualDemand and SpecifiedDemandProfile
    # ### Calculate SpecifiedAnnualDemand and SpecifiedDemandProfile
    
    sp_demand_df = demand_df[
        [x for x in demand_df.columns if x in demand_nodes or x == "TIMESLICE"]
    ]
    
    sp_demand_df = pd.melt(
        sp_demand_df,
        id_vars="TIMESLICE",
        value_vars=demand_nodes,
        var_name="node",
        value_name="demand",
    )
    
    sp_demand_df = sp_demand_df.groupby(["TIMESLICE", "node"], as_index=False).sum()
    
    # Calculate SpecifiedAnnualDemand
    total_demand_df = (
        sp_demand_df.drop(columns="TIMESLICE").groupby("node", as_index=False).sum()
    )
    
    total_demand_df.rename({"demand": "total_demand"}, axis=1, inplace=True)
    
    sp_demand_df = sp_demand_df.join(total_demand_df.set_index("node"), on="node")
    
    # Calculate SpecifiedDemandProfile
    
    sp_demand_df["VALUE"] = sp_demand_df["demand"] / sp_demand_df["total_demand"]
    
    
    # Filter out country aggregate values for countries with multiple nodes
    country_with_nodes = list(
        (sp_demand_df.loc[sp_demand_df["node"].str.len() > 6, "node"].str[:-3].unique())
    )
    
    sp_demand_df = sp_demand_df.loc[~(sp_demand_df["node"].isin(country_with_nodes))]
    
    
    # Rename COMMODITY based on naming convention.
    # Add 'XX' for countries without multiple nodes
    sp_demand_df.loc[sp_demand_df["node"].str.len() == 5, "FUEL"] = (
        "ELC" + sp_demand_df["node"] + "02"
    )
    
    sp_demand_df.loc[sp_demand_df["node"].str.len() == 6, "FUEL"] = (
        "ELC" + sp_demand_df["node"].str.split("-").str[1:].str.join("") + "XX02"
    )
    
    sp_demand_df.loc[sp_demand_df["node"].str.len() > 6, "FUEL"] = (
        "ELC" + sp_demand_df["node"].str.split("-").str[1:].str.join("") + "02"
    )
    
    # In case custom data is provided only keep the custom data
    sp_demand_df.drop_duplicates(subset=['TIMESLICE', 'FUEL'], keep = 'last', inplace = True)
    
    # Create master table for SpecifiedDemandProfile
    sp_demand_df_final = pd.DataFrame(
        list(
            itertools.product(
                sp_demand_df["TIMESLICE"].unique(), sp_demand_df["FUEL"].unique(), years
            )
        ),
        columns=["TIMESLICE", "FUEL", "YEAR"],
    )
    sp_demand_df_final = sp_demand_df_final.join(
        sp_demand_df.set_index(["TIMESLICE", "FUEL"]), on=["TIMESLICE", "FUEL"]
    )
    
    # Add 'REGION' column and fill 'GLOBAL' throughout
    sp_demand_df_final["REGION"] = "GLOBAL"
    
    total_demand_df_final = (
        sp_demand_df_final.groupby(["REGION", "FUEL", "YEAR"], as_index=False)[
            "total_demand"
        ]
        .agg("mean")
        .rename({"total_demand": "VALUE"}, axis=1)
    )
    
    # Convert SpecifiedAnnualDemand to required units
    total_demand_df_final["VALUE"] = total_demand_df_final["VALUE"].mul(3.6 * 1e-6)
    
    # Generate SpecifiedDemandProfile.csv file
    sp_demand_df_final["VALUE"] = sp_demand_df_final["VALUE"].round(2)
    sp_demand_df_final = sp_demand_df_final[
        ["REGION", "FUEL", "TIMESLICE", "YEAR", "VALUE"]
    ].dropna()
    
    # sp_demand_df_final = apply_dtypes(sp_demand_df_final, "SpecifiedDemandProfile")
    sp_demand_df_final.drop_duplicates(
        subset=["REGION", "TIMESLICE", "FUEL", "YEAR"], keep="last", inplace=True
    )
    
    sp_demand_df_final.to_csv(
        os.path.join(output_data_dir, "SpecifiedDemandProfile.csv"), index=None
    )
    
    # CapacityFactor
    
    datetime_ts_df = demand_df[["Datetime", "TIMESLICE"]]
    
    def capacity_factor(df):
        df["Datetime"] = pd.to_datetime(df["Datetime"], format="%d/%m/%Y %H:%M")
        capfac_df = df.set_index("Datetime").join(
            datetime_ts_df.set_index("Datetime"), on="Datetime"
        )
        capfac_nodes = [x for x in capfac_df.columns if x not in ["Datetime", "TIMESLICE"]]
        capfac_df = capfac_df.reset_index().drop("Datetime", axis=1)
        capfac_df = pd.melt(
            capfac_df,
            id_vars="TIMESLICE",
            value_vars=capfac_nodes,
            var_name="node",
            value_name="VALUE",
        )
        capfac_df = capfac_df.groupby(["TIMESLICE", "node"], as_index=False).agg("mean")
        capfac_df["VALUE"] = capfac_df["VALUE"].div(100).round(4)
    
        ## Filter out country aggregate values for countries with multiple nodes
        capfac_df = capfac_df.loc[~(capfac_df["node"].isin(country_with_nodes))]
    
        # Rename COMMODITY based on naming convention.
        # Add 'XX' for countries without multiple nodes
        capfac_df.loc[capfac_df["node"].str.len() == 5, "TECHNOLOGY"] = (
            "PWR" + df.name + capfac_df["node"] + "01"
        )
        
        capfac_df.loc[capfac_df["node"].str.len() == 6, "TECHNOLOGY"] = (
            "PWR" + df.name + capfac_df["node"].str.split("-").str[1:].str.join("") + "XX01"
        )
    
        capfac_df.loc[capfac_df["node"].str.len() > 6, "TECHNOLOGY"] = (
            "PWR" + df.name + capfac_df["node"].str.split("-").str[1:].str.join("") + "01"
        )
        
        # In case custom data is provided only keep the custom data
        capfac_df.drop_duplicates(subset=['TIMESLICE', 'TECHNOLOGY'], keep = 'last', inplace = True)
    
        # Create master table for CapacityFactor
        capfac_df_final = pd.DataFrame(
            list(
                itertools.product(
                    capfac_df["TIMESLICE"].unique(), capfac_df["TECHNOLOGY"].unique(), years
                )
            ),
            columns=["TIMESLICE", "TECHNOLOGY", "YEAR"],
        )
        capfac_df_final = capfac_df_final.join(
            capfac_df.set_index(["TIMESLICE", "TECHNOLOGY"]), on=["TIMESLICE", "TECHNOLOGY"]
        )
    
        # Add 'REGION' column and fill 'GLOBAL' throughout
        capfac_df_final["REGION"] = "GLOBAL"
    
        capfac_df_final = capfac_df_final[
            ["REGION", "TECHNOLOGY", "TIMESLICE", "YEAR", "VALUE"]
        ]
    
        return capfac_df_final
    
    capfacs = []
    for each in [hyd_df_processed, csp_df, spv_df, won_df, wof_df]:
        capfacs.append(capacity_factor(each))
    
    capfac_all_df = pd.concat(capfacs).reset_index(drop=True)
    
    capfac_all_df.drop_duplicates(
        subset=["REGION", "TECHNOLOGY", "TIMESLICE", "YEAR"], keep="last", inplace=True
    )
    capfac_all_df.to_csv(os.path.join(output_data_dir, "CapacityFactor.csv"), index=None)
    
    # Create csv for TIMESLICE
    
    # ## Create csv for TIMESLICE
    time_slice_list = list(demand_df["TIMESLICE"].unique())
    time_slice_df = pd.DataFrame(time_slice_list, columns=["VALUE"]).astype(
        SET_DTYPES["TIMESLICE"]
    )
    time_slice_df.to_csv(os.path.join(output_data_dir, "TIMESLICE.csv"), index=None)
    
    demand_nodes = list(set(list(sp_demand_df_final["FUEL"].str[3:8])))
    
    # Create Conversionls, Conversionld, and Conversionlh
    
    # Conversionls
    df_ls = pd.DataFrame(
        list(itertools.product(time_slice_list, list(range(1, len(seasons) + 1)))),
        columns=["TIMESLICE", "SEASON"],
    )
    df_ls.loc[df_ls["TIMESLICE"].str[1:2].astype(int) == df_ls["SEASON"], "VALUE"] = 1
    df_ls.fillna(0, inplace=True)
    df_ls.to_csv(os.path.join(output_data_dir, "Conversionls.csv"), index=None)
    
    df_season_set = pd.DataFrame(list(range(1, len(seasons) + 1)), columns=["VALUE"])
    df_season_set.to_csv(os.path.join(output_data_dir, "SEASON.csv"), index=None)
    
    # Conversionld
    df_ld = pd.DataFrame(
        list(itertools.product(time_slice_list, [1])), columns=["TIMESLICE", "DAYTYPE"]
    )
    df_ld["VALUE"] = 1
    df_ld.fillna(0, inplace=True)
    df_ld.to_csv(os.path.join(output_data_dir, "Conversionld.csv"), index=None)
    df_daytype_set = pd.DataFrame([1], columns=["VALUE"])
    df_daytype_set.to_csv(os.path.join(output_data_dir, "DAYTYPE.csv"), index=None)
    
    # Conversionlh
    df_lh = pd.DataFrame(
        list(itertools.product(time_slice_list, list(range(1, len(dayparts) + 1)))),
        columns=["TIMESLICE", "DAILYTIMEBRACKET"],
    )
    df_lh.loc[
        df_lh["TIMESLICE"].str[3:].astype(int) == df_lh["DAILYTIMEBRACKET"], "VALUE"
    ] = 1
    df_lh.fillna(0, inplace=True)
    df_lh.to_csv(os.path.join(output_data_dir, "Conversionlh.csv"), index=None)
    df_dayparts_set = pd.DataFrame(list(range(1, len(dayparts) + 1)), columns=["VALUE"])
    df_dayparts_set.to_csv(
        os.path.join(output_data_dir, "DAILYTIMEBRACKET.csv"), index=None
    )
    
    # Daysplit
    
    daysplit = {}
    for dp, hr in dayparts_raw.items():
        daysplit[int(dp[1:])] = (hr[1] - hr[0]) / 8760
    
    df_daysplit = pd.DataFrame(
        itertools.product(list(range(1, len(dayparts) + 1)), years),
        columns=["DAILYTIMEBRACKET", "YEAR"],
    )
    df_daysplit["VALUE"] = df_daysplit["DAILYTIMEBRACKET"].map(daysplit)
    df_daysplit = df_daysplit[["DAILYTIMEBRACKET", "YEAR", "VALUE"]]
    df_daysplit["VALUE"] = df_daysplit["VALUE"].round(4)
    df_daysplit.to_csv(os.path.join(output_data_dir, "DaySplit.csv"), index=None)

if __name__ == "__main__":
    
    if "snakemake" in globals():
        start_year = snakemake.params.start_year
        end_year = snakemake.params.end_year
        region_name = snakemake.params.region_name
        geographic_scope = snakemake.params.geographic_scope
        output_data_dir = snakemake.params.output_data_dir
        input_data_dir = snakemake.params.input_data_dir
        output_dir = snakemake.params.output_dir
        input_dir = snakemake.params.input_dir
        custom_nodes_dir = snakemake.params.input_data_dir
        daytype = snakemake.params.daytype
        seasons = snakemake.params.seasons
        dayparts = snakemake.params.dayparts
        timeshift = snakemake.params.timeshift
        
        plexos_demand = pd.read_csv(snakemake.input.plexos_demand)
        plexos_csp_2015 = pd.read_csv(snakemake.input.plexos_csp_2015)
        plexos_spv_2015 = pd.read_csv(snakemake.input.plexos_spv_2015)
        plexos_hyd_2015 = pd.read_csv(snakemake.input.plexos_hyd_2015)
        plexos_won_2015 = pd.read_csv(snakemake.input.plexos_won_2015)
        plexos_wof_2015 = pd.read_csv(snakemake.input.plexos_wof_2015)
        custom_specified_demand_profiles = pd.read_csv(snakemake.input.custom_specified_demand_profiles)
        custom_csp_profiles = pd.read_csv(snakemake.input.custom_csp_profiles)
        custom_hyd_profiles = pd.read_csv(snakemake.input.custom_hyd_profiles)
        custom_spv_profiles = pd.read_csv(snakemake.input.custom_spv_profiles)
        custom_wof_profiles = pd.read_csv(snakemake.input.custom_wof_profiles)
        custom_won_profiles = pd.read_csv(snakemake.input.custom_won_profiles)
        
    # The below else statement defines variables if the 'powerplant/main' script is to be run locally
    # outside the snakemake workflow. This is relevant for testing purposes only! User inputs when running 
    # the full workflow need to be defined in the config file. 
            
    else:      
        start_year = 2021
        end_year = 2050
        region_name = 'GLOBAL'
        geographic_scope = ['BTN', 'IND']
        output_data_dir = 'results/data'
        input_data_dir = 'resources/data/default'
        output_dir = 'results'
        input_dir = 'resources'
        custom_nodes_dir = 'resources/data/custom'
        daytype = False
        seasons =   {'S1': [1, 2, 3, 4, 5, 6], 
                     'S2': [7, 8, 9, 10, 11, 12]}
        dayparts =   {'D1': [1, 7],
                      'D2': [7, 13],
                      'D3': [13, 19],
                      'D4': [19, 25]}
        timeshift = 0
        
        plexos_demand = pd.read_csv(os.path.join(input_data_dir, 'All_Demand_UTC_2015.csv'), encoding="latin-1")
        plexos_csp_2015 = pd.read_csv(os.path.join(input_data_dir, 'CSP 2015.csv'), encoding="latin-1")
        plexos_spv_2015 = pd.read_csv(os.path.join(input_data_dir, 'SolarPV 2015.csv'), encoding="latin-1")
        plexos_hyd_2015 = pd.read_csv(os.path.join(input_data_dir, 'Hydro_Monthly_Profiles (15 year average).csv'), encoding="latin-1")
        plexos_won_2015 = pd.read_csv(os.path.join(input_data_dir, 'Won 2015.csv'), encoding="latin-1")
        plexos_wof_2015 = pd.read_csv(os.path.join(input_data_dir, 'Woff 2015.csv'), encoding="latin-1")
        custom_specified_demand_profiles = pd.read_csv(os.path.join(custom_nodes_dir, 'specified_demand_profile.csv'))
        custom_csp_profiles = pd.read_csv(os.path.join(custom_nodes_dir, 'RE_profiles_CSP.csv'), encoding="latin-1")
        custom_hyd_profiles = pd.read_csv(os.path.join(custom_nodes_dir, 'RE_profiles_HYD.csv'), encoding="latin-1")
        custom_spv_profiles = pd.read_csv(os.path.join(custom_nodes_dir, 'RE_profiles_SPV.csv'), encoding="latin-1")
        custom_wof_profiles = pd.read_csv(os.path.join(custom_nodes_dir, 'RE_profiles_WOF.csv'), encoding="latin-1")
        custom_won_profiles = pd.read_csv(os.path.join(custom_nodes_dir, 'RE_profiles_WON.csv'), encoding="latin-1")      

    # SET INPUT DATA
    input_data = {
        "demand_df" : plexos_demand,
        "csp_df" : plexos_csp_2015,
        "spv_df" : plexos_spv_2015,
        "hyd_df" : plexos_hyd_2015,
        "won_df" : plexos_won_2015,
        "wof_df" : plexos_wof_2015,
        "custom_sp_demand_profile" : custom_specified_demand_profiles,
        "csp_df_custom" : custom_csp_profiles,
        "hyd_df_custom" : custom_hyd_profiles,
        "spv_df_custom" : custom_spv_profiles,
        "wof_df_custom" : custom_wof_profiles,
        "won_df_custom" : custom_won_profiles,
        "seasons_raw" : seasons,
        "dayparts_raw": dayparts,
    }
    
    # CALL MAIN
    main(**input_data)