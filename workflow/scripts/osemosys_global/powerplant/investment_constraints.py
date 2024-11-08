"""Function to set min and max capacity investment constraints."""

import pandas as pd

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

    df_max_capacity = pd.DataFrame(
        out_data, columns=["REGION", "TECHNOLOGY", "YEAR", "VALUE"]
    )
    df_max_capacity.dropna(inplace=True)
    
    return df_max_capacity

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
    df_max_cap_invest = pd.DataFrame(max_cap_invest_data,
                                    columns = ['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                                    )
    df_max_cap_invest = apply_dtypes(df_max_cap_invest, "TotalAnnualMaxCapacityInvestment")

    df_min_cap_invest = pd.DataFrame(columns = ['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']
                                    )
    df_min_cap_invest = apply_dtypes(df_min_cap_invest, "TotalAnnualMinCapacityInvestment")

    return df_max_cap_invest, df_min_cap_invest