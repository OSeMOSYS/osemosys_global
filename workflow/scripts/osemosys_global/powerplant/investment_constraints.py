"""Functions to set powerplant investment constraints."""
import pandas as pd
import itertools

from data import(
    get_years,
    get_max_value_per_technology
    )

from utils import apply_dtypes

def set_renewable_limits(res_limits, tech_code_dict,
                         custom_nodes, custom_nodes_res_limits,
                         residual_capacity, start_year, 
                         end_year, region_name):
    
    years = get_years(start_year, end_year)

    # GET CAPACITY LIMITS FROM PLEXOS WORLD
    # This is capacity ADDITION limits, not total max capacity limits

    df_reslimit_units = res_limits.loc[
        (res_limits["child_object"].str.contains("|".join(tech_code_dict.keys())))
        & (res_limits["property"] == "Max Units Built")
        & (res_limits["scenario"].str.contains("Base"))
        & (res_limits["child_class"] == "Generator")
    ].set_index("child_object")

    df_reslimit_capacity = res_limits.loc[
        (res_limits["child_object"].str.contains("|".join(tech_code_dict.keys())))
        & (res_limits["property"] == "Max Capacity")
        & (res_limits["child_class"] == "Generator")
    ].set_index("child_object")

    df_reslimit_final = pd.DataFrame(
        df_reslimit_capacity["value"] * df_reslimit_units["value"] / 1000
    ).rename(columns={"value": "VALUE"})
    df_reslimit_final["node"], df_reslimit_final["powerplant"] = (
        df_reslimit_final.index.str.rsplit(pat="|",n=1).str[1],
        df_reslimit_final.index.str.rsplit(pat="|",n=1).str[0],
    )
    df_reslimit_final["powerplant"] = df_reslimit_final["powerplant"].map(tech_code_dict)

    df_reslimit_final.loc[df_reslimit_final["node"].str.len() <= 6, "node_code"] = (
        df_reslimit_final["node"].str.split("-").str[1:].str.join("") + "XX"
    )
    df_reslimit_final.loc[df_reslimit_final["node"].str.len() > 6, "node_code"] = (
        df_reslimit_final["node"].str.split("-").str[1:].str.join("")
    )

    df_reslimit_final["TECHNOLOGY"] = (
        "PWR" + df_reslimit_final["powerplant"] + df_reslimit_final["node_code"] + "01"
    )
    cap_addition_limit = df_reslimit_final.set_index("TECHNOLOGY").to_dict()["VALUE"]

    if custom_nodes:

        custom_nodes_res_limits["TECHNOLOGY"] = (
            "PWR"
            + custom_nodes_res_limits["FUEL_TYPE"]
            + custom_nodes_res_limits["CUSTOM_NODE"]
            + "01"
        )
        re_potentials_dict = dict(
            zip(
                list(custom_nodes_res_limits["TECHNOLOGY"]),
                list(custom_nodes_res_limits["CAPACITY"]),
            )
        )
        cap_addition_limit.update(re_potentials_dict)

    # GET RESIDUAL CAPACITY VALUES

   
    residual_capacity["VALUE"] = residual_capacity.loc[:, "VALUE"].round(4)
    df_res_cap = residual_capacity.loc[
        residual_capacity["TECHNOLOGY"].str[3:6].isin(list(tech_code_dict.values()))
    ]
    df_res_cap = get_max_value_per_technology(df_res_cap)
    res_cap = df_res_cap.set_index("TECHNOLOGY").to_dict()["VALUE"]

    # CALCULATE AND FORMAT DATA

    out_data = []
    for tech, maxad in cap_addition_limit.items():
        try:
            max_capacity = float(res_cap[tech]) + float(maxad)
        except KeyError:
            max_capacity = float(maxad)

        # Add 0.0002 to enusre there is no rounding mismathch between total
        # annual max capacity and residual capacity
        max_capacity = round(max_capacity, 4) + 0.0002

        for year in years:
            out_data.append([region_name, tech, year, max_capacity])

    df_max_cap_investacity = pd.DataFrame(
        out_data, columns=["REGION", "TECHNOLOGY", "YEAR", "VALUE"]
    )
    df_max_cap_investacity.dropna(inplace=True)
    
    return df_max_cap_investacity

def cap_investment_constraints(df_iar_final, no_investment_techs,
                               start_year, end_year, region_name):

    # Create totalAnnualMaxCapacityInvestment data 

    # Do not allow capacity investment for all PWRxxxxxxxx00 technologies.
    max_cap_invest_techs = list(set(
        df_iar_final.loc[df_iar_final['TECHNOLOGY'].str.endswith('00')]['TECHNOLOGY'].tolist()))
    max_cap_invest_data = []
    for tech in max_cap_invest_techs:
        for year in get_years(start_year, end_year): 
            max_cap_invest_data.append([region_name, tech, year, 0])

    # Do not allow investment for all xxxABCxxxxxxx technologies
    
    if not no_investment_techs:
        no_investment_techs = [] # Change from None type to empty list
    max_cap_invest_techs = list(set(df_iar_final.loc[
        df_iar_final['TECHNOLOGY'].str[3:6].isin(no_investment_techs)][
        'TECHNOLOGY'].tolist()))
    for tech in max_cap_invest_techs:
        for year in get_years(start_year, end_year): 
            max_cap_invest_data.append([region_name, tech, year, 0])
    
    # Save totalAnnualMaxCapacityInvestment
    df_max_cap_invest_invest = pd.DataFrame(max_cap_invest_data,
                                    columns = ['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                                    )
    df_max_cap_invest_invest = apply_dtypes(df_max_cap_invest_invest, "TotalAnnualMaxCapacityInvestment")

    df_min_cap_invest = pd.DataFrame(columns = ['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                                    )
    df_min_cap_invest = apply_dtypes(df_min_cap_invest, "TotalAnnualMinCapacityInvestment")

    return df_max_cap_invest_invest, df_min_cap_invest

def set_build_rates(build_rates, tech_set, max_cap_invest, 
                    max_cap, start_year, end_year, region_name):
    
    years = get_years(start_year, end_year)
    
    max_build_df = build_rates.loc[
        build_rates.index.repeat(
            build_rates["END_YEAR"] + 1 - build_rates["START_YEAR"])]
    
    max_build_df["YEAR"] = (
        max_build_df.groupby(level=0).cumcount() + max_build_df["START_YEAR"]
    )
    
    max_build_df = max_build_df.reset_index(drop=True)
    max_build_df = max_build_df[["TYPE", "METHOD", "MAX_BUILD", "YEAR", "COUNTRY"]]
    max_build_df["TYPE"] = max_build_df["TYPE"].str[0:3]
    # Create a list of powerplant technologies
    pwr_tech_list = [x for x in list(tech_set["VALUE"]) if x.startswith("PWR")]

    # Create scaffold dataframe with all powerplant technologies for all years
    df_techs = pd.DataFrame(
        list(itertools.product(pwr_tech_list, years)), columns=["TECHNOLOGY", "YEAR"]
    )

    # Filter out technologies for which a max. capacity investment has already
    # been set
    max_cap_inv_techs = list(max_cap_invest["TECHNOLOGY"].unique())
    df_techs = df_techs[~(df_techs["TECHNOLOGY"].isin(max_cap_inv_techs))]

    # Create dataframe of max capacity by technology
    df_max_cap_invest = max_cap.loc[max_cap["TECHNOLOGY"].str.startswith("PWR")]

    df_max_cap_invest = pd.merge(
        left=df_techs, right=df_max_cap_invest, on=["TECHNOLOGY", "YEAR"], how="left"
    )
    df_max_cap_invest["REGION"] = region_name
    df_max_cap_invest["TYPE"] = df_max_cap_invest["TECHNOLOGY"].str[3:6]
    df_max_cap_invest["COUNTRY"] = df_max_cap_invest["TECHNOLOGY"].str[6:9]
    df_max_cap_invest = pd.merge(
        left=df_max_cap_invest, right=max_build_df, on=["TYPE", "COUNTRY", "YEAR"], how="left"
    )
    df_max_cap_invest.loc[df_max_cap_invest["METHOD"].isin(["ABS"]), "VALUE"] = df_max_cap_invest[
        "MAX_BUILD"
    ]
    df_max_cap_invest.dropna(inplace=True)
    df_max_cap_invest.loc[df_max_cap_invest["METHOD"].isin(["PCT"]), "VALUE"] = (
        df_max_cap_invest["VALUE"] * df_max_cap_invest["MAX_BUILD"] / 100
    )

    df_max_cap_invest = df_max_cap_invest[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    df_max_cap_invest = pd.concat([max_cap_invest, df_max_cap_invest], ignore_index=True)
    df_max_cap_invest["VALUE"] = df_max_cap_invest["VALUE"].astype(float).round(3)
    
    return df_max_cap_invest

def set_fossil_capacity_constraints(tech_list, fossil_targets, 
                                    df_max_capacity, df_min_capacity,
                                    df_residual_capacity, region_name):
    
    df_min_capacity.columns = ["REGION", "TECHNOLOGY", "YEAR", "VALUE"]
    
    if fossil_targets:
        
        # Loop through targets and apply
        for param in fossil_targets:
            region = param[0]
            fuel = param[1]
            technology = 'PWR' + fuel + region
            first_year = param[2]
            last_year = param[3]
            sense = param[4]
            value = param[5]
            
            years = get_years(first_year, last_year)
            
            # Pull residual capacity for OCG/CCG '00' technologies
            res_cap = df_residual_capacity.loc[
                (df_residual_capacity['TECHNOLOGY'] == technology + '00')
                ]
            
            res_cap.loc[
                res_cap['TECHNOLOGY'].str.contains('00'), 'TECHNOLOGY'
                ] = res_cap['TECHNOLOGY'].str.replace('00', '01') 
                                  
            if sense == 'ABS' or sense == 'MAX':
                    
                    custom_max_cap = pd.DataFrame(
                        list(itertools.product([region_name], [technology + '01'], 
                                               years, [value])), 
                        columns=["REGION", "TECHNOLOGY", "YEAR", "VALUE"]
                    )
                    
                    # Substract residual capacity for OCG/CCG '00' technologies if applicable
                    if not res_cap.empty:
                        custom_max_cap = pd.merge(custom_max_cap, res_cap[['YEAR', 'VALUE']], 
                                                  left_on = 'YEAR', right_on = 'YEAR')
                        custom_max_cap['VALUE'] = custom_max_cap['VALUE_x'
                                                                 ] - custom_max_cap['VALUE_y']
                    
                    df_max_capacity = pd.concat([df_max_capacity, custom_max_cap], 
                                                join = 'inner').reset_index(drop = True)
                          
            if sense == 'ABS' or sense == 'MIN':
                    
                    custom_min_cap = pd.DataFrame(
                        list(itertools.product([region_name], [technology + '01'], 
                                               years, [value])), 
                        columns=["REGION", "TECHNOLOGY", "YEAR", "VALUE"]
                    )
                    
                    # Substract residual capacity for OCG/CCG '00' technologies if applicable
                    if not res_cap.empty:
                        custom_min_cap = pd.merge(custom_min_cap, res_cap[['YEAR', 'VALUE']], 
                                                  left_on = 'YEAR', right_on = 'YEAR')
                        custom_min_cap['VALUE'] = custom_min_cap['VALUE_x'
                                                                 ] - custom_min_cap['VALUE_y']
                    
                    df_min_capacity = pd.concat([df_min_capacity, custom_min_cap],
                                                join = 'inner').reset_index(drop = True)
    
    return df_max_capacity, df_min_capacity