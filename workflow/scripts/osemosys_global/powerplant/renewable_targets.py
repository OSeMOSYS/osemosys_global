"""Function to set operational life."""
import pandas as pd
import itertools

from data import get_years

def apply_re_pct_targets(re_targets, remove_nodes, oar_df, renewables_list, 
                     fuel_set, annual_demand_df, region_name):
    
    """Apply Renewable Energy targets by country and year"""
    
    if not remove_nodes:
        remove_nodes = []
        
    country_dict = {}
    target_techs_dict = {}
    type_dict = {}
    first_year_dict = {}
    final_year_dict = {}
    value_dict = {}

    # Get transmission parameters from config file.
    for target, target_params in re_targets.items():
        country_dict[target] = target_params[0]
        target_techs_dict[target] = target_params[1]
        type_dict[target] = target_params[2]
        first_year_dict[target] = target_params[3]
        final_year_dict[target] = target_params[4]
        value_dict[target] = target_params[5]
    
    accumulated_annual_demand = pd.DataFrame(columns = ["REGION", "FUEL", "YEAR", "VALUE"])
    
    if not re_targets is None:
        for target in re_targets.keys():
            if type_dict[target] == 'PCT':
                if not target_techs_dict[target]:
                    tech_filter = renewables_list
                else:
                    tech_filter = target_techs_dict[target]

                re_df = oar_df.loc[
                    (oar_df["TECHNOLOGY"].str.startswith("PWR"))
                    & (oar_df["TECHNOLOGY"].str[3:6].isin(tech_filter))
                    & (oar_df["TECHNOLOGY"].str[6:9] == country_dict[target])
                ].copy()
                
                # Create dummy commodity starting with 'REN' for renewables
                re_df["FUEL"] = target + country_dict[target]
                oar_df = pd.concat([oar_df, re_df]).drop_duplicates()
            
                # Create list of fuels
                fuels_ren_df = re_df[["FUEL"]].copy()
                fuels_ren_df.rename(columns={"FUEL": "VALUE"}, inplace=True)
    
                fuel_set = pd.concat([fuel_set, fuels_ren_df]).dropna()
                fuel_set.drop_duplicates(inplace=True)
     
                # Create dataframe template to calculate SpecifiedAnnualDemand
                target_years = get_years(first_year_dict[target], final_year_dict[target])
                re_targets_df = pd.DataFrame(
                    list(itertools.product(target + country_dict[target], target_years)), 
                    columns=["FUEL", "YEAR"]
                )
                re_targets_df["COUNTRY"] = country_dict[target]
                re_targets_df = re_targets_df[["COUNTRY", "YEAR"]]
    
                re_targets_df['VALUE'] = value_dict[target] / 100
                re_targets_df.drop_duplicates(inplace=True)
    
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
                sp_demand_df = annual_demand_df.loc[
                    ~(annual_demand_df["FUEL"].str[3:8].isin(remove_nodes))
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
                re_targets_df["FUEL"] = target + re_targets_df["COUNTRY"]
                re_targets_df["REGION"] = region_name
                re_targets_df = re_targets_df[["REGION", "FUEL", "YEAR", "VALUE"]]
                
                if accumulated_annual_demand.empty:
                
                    accumulated_annual_demand = re_targets_df
                
                else:
                    accumulated_annual_demand = pd.concat([accumulated_annual_demand, re_targets_df])

    return fuel_set, oar_df, accumulated_annual_demand

def apply_re_abs_targets(re_targets, remove_nodes, oar_df, renewables_list, 
                     fuel_set, annual_demand_df, region_name):
    
    a = ''