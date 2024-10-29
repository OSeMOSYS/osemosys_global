"""Function to calculate transmission costs."""
import pandas as pd

from data import(
    get_years,
    set_transmission_tech_groups,
    calculate_transmission_distances
    )

def get_transmission_costs(df_exist_corrected, df_planned_corrected,
                           df_pp_capex_base, df_pp_fix_base,
                           df_pp_var_base, centerpoints_dict, 
                           trn_param, start_year, end_year, 
                           region_name, subsea_lines):
    '''Gets electrical transmission capital, fixed and variable cost per technology. 

    Both the capital costs and fixed cost are written out to avoid having 
    to read in the excel file twice 
    
    Returns: 
        df_capex: DataFrame with the columns 'TECHNOLOGY' and 'VALUE' 
            representing CAPITAL cost in millions of dollars per year. 
        df_fix: DataFrame with the columns 'TECHNOLOGY' and 'VALUE' 
            representing FIXED cost in millions of dollars per year.
        df_var: DataFrame with the columns 'TECHNOLOGY' and 'VALUE' 
            representing VARIABLE cost in millions of dollars/PJ.             
    '''
    line_capex_dict = {}
    converter_capex_dict = {}
    fom_dict = {}
    var_dict = {}    

    # Get transmission parameters from config file.
    for tech, tech_params in trn_param.items():
        line_capex_dict[tech] = tech_params[0]
        converter_capex_dict[tech] = tech_params[1]
        fom_dict[tech] = tech_params[4]
        var_dict[tech] = tech_params[5]        
    
    # Assign technology groups to individual technologies.
    tech_group_dict = set_transmission_tech_groups(df_exist_corrected, 
                                                   df_planned_corrected, 
                                                   centerpoints_dict, 
                                                   trn_param, 
                                                   subsea_lines)
    
    # Set base df including calculated transmission distances per pathway.
    df_capex = calculate_transmission_distances(df_exist_corrected, df_planned_corrected, 
                                          centerpoints_dict)[['TECHNOLOGY', 'distance']]
    
    # Calculate technology and line specific CAPEX costs ($000).
    df_capex['tech_group'] = df_capex['TECHNOLOGY'].map(tech_group_dict)
    
    for group in trn_param.keys():
        df_capex.loc[df_capex['tech_group'] == group, 
                     'VALUE'] = (converter_capex_dict[group] + (
                         line_capex_dict[group] * df_capex['distance'])) / 100

    # Add missing columns to df.
    df_capex['REGION'] = region_name
    df_capex['YEAR'] = [get_years(start_year, end_year)] * len(df_capex)
    df_capex = df_capex.explode('YEAR').reset_index()

    # Set Fixed O&M costs as a function of total CAPEX and set default var costs.
    df_fix = df_capex.copy()
    df_var = df_capex.copy()
    
    for group in trn_param.keys():
        df_fix.loc[df_fix['tech_group'] == group, 
                     'VALUE'] = df_fix['VALUE'] * (fom_dict[group] / 100)
        
        # Variable cost in $/MWh converted to $Million/PJ
        df_var.loc[df_var['tech_group'] == group, 
                   'VALUE'] = round(var_dict[group] / 0.0000036 / 1000000 , 4)
                
        
    df_var['MODE_OF_OPERATION'] = [[1,2] for x in range(len(df_var))]
    df_var = df_var.explode('MODE_OF_OPERATION')
    
    # Filter out required columns.
    df_capex = df_capex[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]
    df_fix = df_fix[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]
    df_var = df_var[['REGION', 'TECHNOLOGY', 'MODE_OF_OPERATION', 'YEAR', 'VALUE']]    
    
    # Add transmission costs dfs to powerplant costs dfs.
    
    df_capex = pd.concat([df_pp_capex_base, df_capex])
    df_fix = pd.concat([df_pp_fix_base, df_fix])
    df_var = pd.concat([df_pp_var_base, df_var])    
    
    return df_capex, df_fix, df_var