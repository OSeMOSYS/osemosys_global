"""Function to set renewable targets."""
import pandas as pd
import itertools

from data import get_years

def apply_re_pct_targets(re_targets, geographic_scope, remove_nodes, oar_df, 
                         renewables_list, fuel_set, specified_demand_df, region_name):
    
    """Apply Renewable Energy targets by country, year and technology in 
    relative terms compared to overall generation (= specified demand)."""
    
    if not remove_nodes:
        remove_nodes = []
        
    node_dict = {}
    target_techs_dict = {}
    type_dict = {}
    first_year_dict = {}
    final_year_dict = {}
    value_dict = {}    
    accumulated_annual_demand = pd.DataFrame(columns = ["REGION", "FUEL", 
                                                        "YEAR", "VALUE"])
    
    # Check if any targets are defined.
    if not re_targets is None:
        # Get transmission parameters from config file.
        for target, target_params in re_targets.items():
            node_dict[target] = target_params[0]
            target_techs_dict[target] = target_params[1]
            type_dict[target] = target_params[2]
            first_year_dict[target] = target_params[3]
            final_year_dict[target] = target_params[4]
            value_dict[target] = target_params[5]
        
        # Per entry check if they are relative targets ('PCT').
        for target in re_targets:
            if type_dict[target] == 'PCT':
                # Check for a technology subset for the target.
                if not target_techs_dict[target]:
                    tech_filter = renewables_list
                else:
                    tech_filter = target_techs_dict[target]
                 
                # Filter OAR df to only keep the relevant country techs.
                if len(str(node_dict[target])) == 3:
                    re_df = oar_df.loc[
                        (oar_df["TECHNOLOGY"].str.startswith("PWR"))
                        & (oar_df["TECHNOLOGY"].str[3:6].isin(tech_filter))
                        & (oar_df["TECHNOLOGY"].str[6:9] == node_dict[target])
                    ].copy()
                
                elif len(str(node_dict[target])) == 5:
                    re_df = oar_df.loc[
                        (oar_df["TECHNOLOGY"].str.startswith("PWR"))
                        & (oar_df["TECHNOLOGY"].str[3:6].isin(tech_filter))
                        & (oar_df["TECHNOLOGY"].str[6:11] == node_dict[target])
                    ].copy()
                    
                elif len(str(node_dict[target])) == 0:
                    re_df = oar_df.loc[
                        (oar_df["TECHNOLOGY"].str.startswith("PWR"))
                        & (oar_df["TECHNOLOGY"].str[3:6].isin(tech_filter))
                        & (oar_df["TECHNOLOGY"].str[6:9].isin(geographic_scope))
                    ].copy()
                
                # Create dummy commodity for the target.
                re_df["FUEL"] = target + node_dict[target]
                
                # Add dummy commodity to OAR.
                oar_df = pd.concat([oar_df, re_df]).drop_duplicates()
            
                # Add dummy commodity to fuel list.
                fuels_ren_df = re_df[["FUEL"]].copy()
                fuels_ren_df.rename(columns={"FUEL": "VALUE"}, inplace=True)
    
                fuel_set = pd.concat([fuel_set, fuels_ren_df]).dropna()
                fuel_set.drop_duplicates(inplace=True)
     
                # Create dataframe template to calculate AccumulatedAnnualDemand
                target_years = get_years(first_year_dict[target], final_year_dict[target])
                re_targets_df = pd.DataFrame(
                    list(itertools.product(target + node_dict[target], target_years)), 
                    columns=["FUEL", "YEAR"]
                )
                re_targets_df["GEO"] = node_dict[target]
                re_targets_df = re_targets_df[["GEO", "YEAR"]]
                
                re_targets_df['VALUE'] = value_dict[target] / 100
                re_targets_df.drop_duplicates(inplace=True)
    
                re_targets_df = re_targets_df.pivot(
                    index=["YEAR"], columns=["GEO"], values="VALUE"
                ).reset_index()
                re_targets_df = re_targets_df.interpolate()
        
                # Drop all columns with only NaN
                re_targets_df.dropna(axis=1, how="all", inplace=True)
        
                # Melt to get 'COUNTRY' column and remove all rows with NaN
                re_targets_df = pd.melt(
                    re_targets_df,
                    id_vars=["YEAR"],
                    value_vars=[x for x in re_targets_df.columns if x not in ["YEAR"]],
                    var_name="GEO",
                    value_name="VALUE",
                )
                re_targets_df.dropna(axis=0, inplace=True)
        
                # Read 'SpecifiedAnnualDemand' as output from the demand_projection rule.
                sp_demand_df = specified_demand_df.loc[
                    ~(specified_demand_df["FUEL"].str[3:8].isin(remove_nodes))
                ]
                  
                if len(str(node_dict[target])) == 3:
                    sp_demand_df["GEO"] = sp_demand_df["FUEL"].str[3:6]
                    
                elif len(str(node_dict[target])) == 5:
                    sp_demand_df["GEO"] = sp_demand_df["FUEL"].str[3:8]
                    
                elif len(str(node_dict[target])) == 0:
                    sp_demand_df["GEO"] = ''
                    sp_demand_df = sp_demand_df.loc[sp_demand_df["FUEL"
                                                                 ].str[3:6].isin(
                                                                     geographic_scope)]
                    
                sp_demand_df = sp_demand_df.groupby(["YEAR", "GEO"], as_index=False)[
                    "VALUE"
                ].sum()
                sp_demand_df.rename(columns={"VALUE": "DEMAND"}, inplace=True)
                
                # Merge RE targets and SpecifiedAnnualDemand tables
                re_targets_df = pd.merge(
                    re_targets_df, sp_demand_df, how="left", on=["GEO", "YEAR"]
                )
                
                # Calculate the target for the dummy RES commodity
                re_targets_df["VALUE"] = re_targets_df["VALUE"] * re_targets_df["DEMAND"]
                re_targets_df["FUEL"] = target + re_targets_df["GEO"]
                re_targets_df["REGION"] = region_name
                re_targets_df = re_targets_df[["REGION", "FUEL", "YEAR", "VALUE"]]
                
                if accumulated_annual_demand.empty:
                
                    accumulated_annual_demand = re_targets_df
                
                else:
                    accumulated_annual_demand = pd.concat([
                        accumulated_annual_demand, 
                        re_targets_df])

    return fuel_set, oar_df, accumulated_annual_demand

def apply_re_abs_targets(re_targets, remove_nodes, region_name):
    
    """Apply Renewable Energy targets by country, year and technology
    in absolute terms (GW)."""
    
    if not remove_nodes:
        remove_nodes = []
        
    node_dict = {}
    target_techs_dict = {}
    type_dict = {}
    first_year_dict = {}
    final_year_dict = {}
    value_dict = {}
        
    total_annual_min_capacity = pd.DataFrame(columns = ["REGION", "TECHNOLOGY", 
                                                        "YEAR", "VALUE"])
    # Check if any targets are defined.
    if not re_targets is None:
        # Get transmission parameters from config file.
        for target, target_params in re_targets.items():
            node_dict[target] = target_params[0]
            target_techs_dict[target] = target_params[1]
            type_dict[target] = target_params[2]
            first_year_dict[target] = target_params[3]
            final_year_dict[target] = target_params[4]
            value_dict[target] = target_params[5]

        # Per entry check if they are absolute targets ('ABS').
        for target in re_targets:
            if type_dict[target] == 'ABS':
                target_years = get_years(first_year_dict[target], 
                                         final_year_dict[target])
                # Check to make sure only one tech is defined per 'ABS' target.
                if len(target_techs_dict[target]) != 1:
                    raise Exception(
                        f'Only one tech can be defined per ABS target as entered \
                            in the re_targets config parameter. {target} is \
                                incorrectly set.')

                # Check to make sure 'ABS' target is defined at the nodal level.
                if len(str(node_dict[target])) != 5:
                    raise Exception(
                        f'ABS targets can only be set at the nodal level \
                            in the re_targets config parameter (e.g. INDSO). \
                                {target} is incorrectly set.')
                
                # Set target specific data
                data = pd.DataFrame(
                    list(itertools.product(target_years)), 
                    columns=[["YEAR"]])
                
                data['REGION'] = region_name
                data['TECHNOLOGY'] = 'PWR' + target_techs_dict[
                    target][0] + node_dict[target] + '01'
                data['VALUE'] = value_dict[target]
                
                # Combine data
                if total_annual_min_capacity.empty:
                
                    total_annual_min_capacity = data
                
                else:
                    total_annual_min_capacity = pd.concat([
                        total_annual_min_capacity, 
                        data])
    
    total_annual_min_capacity = total_annual_min_capacity[["REGION", "TECHNOLOGY",
                                                          "YEAR", "VALUE"]]
    
    return total_annual_min_capacity  