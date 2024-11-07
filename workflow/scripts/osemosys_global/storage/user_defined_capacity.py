"""Function to integrate user defined capacities for storage."""

import pandas as pd
import itertools

from data import get_years

def set_user_defined_capacity_sto(tech_capacity_sto,
                                  storage_param,
                                  op_life_dict, 
                                  min_cap_invest_base, 
                                  max_cap_invest_base, 
                                  res_cap_base,
                                  res_cap_sto_base,
                                  df_oar_base,
                                  cap_cost_base,
                                  cap_cost_sto_base, 
                                  fix_cost_base, 
                                  var_cost_base,
                                  start_year, 
                                  end_year, 
                                  region_name
                                  ):
    
    techCapacity_sto = []
    first_year_expansion_dict = {}
    build_rate_dict = {}
    capex_dict = {}
    fom_dict = {}
    var_dict = {}    
    build_year_dict = {}
    efficiency_dict = {}
    duration_dict = {}    

    for idx, tech_params in tech_capacity_sto.items():
        techCapacity_sto.append([idx, tech_params[0], tech_params[1], tech_params[2]])
        build_year_dict[idx] = tech_params[2]
        first_year_expansion_dict[idx] = tech_params[3]
        build_rate_dict[idx] = tech_params[4]
        capex_dict[idx] = tech_params[5]
        fom_dict[idx] = tech_params[6]
        var_dict[idx] = tech_params[7]        
        efficiency_dict[idx] = tech_params[8]
        
    # Set baseline duration as defined in the config file.
    for tech, tech_params in storage_param.items():
        duration_dict[tech] = tech_params[4]
                
    tech_capacity_sto_df = pd.DataFrame(techCapacity_sto,
                                    columns=['idx', 'TECHNOLOGY', 'VALUE', 'YEAR'])
    
    tech_capacity_sto_df['REGION'] = region_name
    
    df_min_cap_inv = pd.concat([min_cap_invest_base, tech_capacity_sto_df])
    df_min_cap_inv.drop_duplicates(inplace=True)
    
    max_cap_techs_df = pd.DataFrame(list(itertools.product(list(tech_capacity_sto_df['idx']),
                                             get_years(start_year, end_year))
                           ),
                      columns = ['idx', 
                                 'YEAR']
                      )
    
    max_cap_techs_df['REGION'], max_cap_techs_df['VALUE'] = region_name, ''
    
    max_cap_techs_df = pd.merge(max_cap_techs_df, df_min_cap_inv[['TECHNOLOGY', 'idx']],
                  how='left',
                  on=['idx'])
    
    max_cap_techs_df = max_cap_techs_df[['idx', 'REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]

    for idx, tech_params in tech_capacity_sto.items():
        max_cap_techs_df.loc[max_cap_techs_df['idx'] == idx, 'BUILD_YEAR'] = build_year_dict.get(idx)
        max_cap_techs_df.loc[max_cap_techs_df['idx'] == idx, 'FIRST_YEAR'] = first_year_expansion_dict.get(idx)
        max_cap_techs_df.loc[max_cap_techs_df['idx'] == idx, 'MAX_BUILD'] = build_rate_dict.get(idx)
            
    max_cap_techs_df.loc[(max_cap_techs_df['YEAR']>=max_cap_techs_df['FIRST_YEAR']),
           'VALUE'] = max_cap_techs_df['MAX_BUILD']
        
    with pd.option_context("future.no_silent_downcasting", True):
        max_cap_techs_df['VALUE'] = max_cap_techs_df[
            'VALUE'].replace(r'^\s*$', 0, regex = True).infer_objects(copy=False)

    max_cap_techs_df = max_cap_techs_df[['REGION',
                                         'TECHNOLOGY',
                                         'YEAR',
                                         'VALUE']]
    
    # Append existing TotalAnnualMaxCapacityInvestment data with MAX_BUILD
    df_max_cap_inv = pd.concat([max_cap_invest_base, max_cap_techs_df]).drop_duplicates()

    # Print TotalAnnualMaxCapacityInvestment.csv with MAX_BUILD
    df_max_cap_inv.drop_duplicates(subset=['REGION', 
                                           'TECHNOLOGY',
                                           'YEAR'],
                                   keep='last',
                                   inplace=True)
        
    # For technologies with start year before model start year, add to ResidualCapacity
    df_res_cap_ud = df_min_cap_inv.copy().loc[df_min_cap_inv.copy()['YEAR'] < min(get_years(start_year, end_year))]
    df_res_cap_ud.rename(columns={'YEAR':'START_YEAR'},
                         inplace=True)
    df_res_cap_ud_final = pd.DataFrame(list(itertools.product(df_res_cap_ud['TECHNOLOGY'].unique(),
                                                              get_years(start_year, end_year))
                                            ),
                                       columns = ['TECHNOLOGY',
                                                  'YEAR']
                                       )
    df_res_cap_ud_final = pd.merge(df_res_cap_ud_final,
                                   df_res_cap_ud,
                                   how='left',
                                   on=['TECHNOLOGY'])
    
    df_res_cap_ud_final['TECH'] = df_res_cap_ud_final['TECHNOLOGY'].str[3:6]
    df_res_cap_ud_final['OP_LIFE'] = df_res_cap_ud_final['TECH'].map(op_life_dict)
    df_res_cap_ud_final['END_YEAR'] = (df_res_cap_ud_final['OP_LIFE'] 
                                       + df_res_cap_ud_final['START_YEAR'])
    df_res_cap_ud_final = df_res_cap_ud_final.loc[df_res_cap_ud_final['YEAR'] 
                                                  >= df_res_cap_ud_final['START_YEAR']]
    df_res_cap_ud_final = df_res_cap_ud_final.loc[df_res_cap_ud_final['YEAR'] 
                                                  <= df_res_cap_ud_final['END_YEAR']]
    df_res_cap_ud_final['REGION'] = region_name
    df_res_cap_ud_final = df_res_cap_ud_final[['REGION',
                                               'TECHNOLOGY',
                                               'YEAR',
                                               'VALUE']]
    
    df_res_cap = pd.concat([res_cap_base, df_res_cap_ud_final 
                            if not df_res_cap_ud_final.empty else None])

    # Group residual capacities in case user defined technology entries already exist.
    df_res_cap = df_res_cap.groupby(['REGION', 'TECHNOLOGY', 'YEAR']
                                    , as_index = False).sum()
    
    # Do the same for ResidualStorageCapacity
    df_res_cap_sto_ud_final = df_res_cap_ud_final.copy().rename(
        columns = {'TECHNOLOGY' : 'STORAGE'}).astype({'VALUE' : float})
    
    df_res_cap_sto_ud_final['STORAGE'] = df_res_cap_sto_ud_final[
        'STORAGE'].str.replace('PWR', '', regex=True)
    
    for tech, duration in duration_dict.items():
        df_res_cap_sto_ud_final.loc[df_res_cap_sto_ud_final[
            'STORAGE'].str.startswith(tech), 'VALUE'] = df_res_cap_sto_ud_final.loc[
                df_res_cap_sto_ud_final['STORAGE'].str.startswith(tech), 'VALUE'] * \
                    duration / 277.777778
                    
    df_res_sto_cap = pd.concat([res_cap_sto_base, df_res_cap_sto_ud_final 
                            if not df_res_cap_sto_ud_final.empty else None])

    # Group residual capacities in case user defined technology entries already exist.
    df_res_sto_cap = df_res_sto_cap.groupby(['REGION', 'STORAGE', 'YEAR']
                                    , as_index = False).sum()                    

    # For technologies with start year at or after model start year, add to 
    # TotalAnnualMinCapacityInvestment      
    df_min_cap_inv = df_min_cap_inv[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']].loc[
        df_min_cap_inv['YEAR'] >= min(get_years(start_year, end_year))]
    
    df_min_cap_inv.drop_duplicates(subset=['REGION', 
                                           'TECHNOLOGY',
                                           'YEAR'],
                                   keep='last',
                                   inplace=True)

    # Update OAR with user-defined efficiencies by storage technology
    df_oar = df_oar_base.copy()

    for idx, tech_params in tech_capacity_sto.items():
        df_oar.loc[df_oar['TECHNOLOGY'] == tech_params[0],
                   'VALUE'] = round(1 / (efficiency_dict[idx] / 100), 3)
    
    # Update CapitalCostStorage with user-defined capex costs by storage technology
    df_cap_cost_sto = cap_cost_sto_base.copy()
    
    ''' Sets capital cost by taking the defined capital cost (m$/GW) divided by the storage 
    duration (=Storage Capacity (GWh)/Storage Power Rating (GW)) to get to GWh values followed 
    by the conversion to PJ (1 GWh = 0.0036 PJ). The resulting costs are divided by 2 to split
    the costs of the overall technology between the charging/discharging component ('PWR') and
    the storage component.'''
    for idx, tech_params in tech_capacity_sto.items():
          df_cap_cost_sto.loc[df_cap_cost_sto['STORAGE'] == 
                         tech_params[0].replace('PWR', ''),
                         'VALUE'] = capex_dict[idx] / duration_dict[
                             tech_params[0][3:6]] / 0.0036 / 2
    
    # Update CapitalCost with user-defined capex costs by storage technology
    df_cap_cost = cap_cost_base.copy()
    
    ''' Sets capital cost by taking the defined capital cost (m$/GW) divided by 2 to split
    the costs of the overall technology between the charging/discharging component ('PWR') and
    the storage component.'''
    for idx, tech_params in tech_capacity_sto.items():
          df_cap_cost.loc[df_cap_cost['TECHNOLOGY'] == tech_params[0],
                         'VALUE'] = capex_dict[idx] / 2
        
    # Update FixedCost with user-defined fixed costs by storage technology
    df_fix_cost = fix_cost_base.copy()

    ''' Sets fixed cost by taking the defined fixed cost (m$/GW/yr) divided by the storage 
    duration (=Storage Capacity (GWh)/Storage Power Rating (GW)) to mimic fixed costs for storage
    as capacities are set in GWh/PJ terms.'''
    for idx, tech_params in tech_capacity_sto.items():
        df_fix_cost.loc[df_fix_cost['TECHNOLOGY'] == tech_params[0],
                   'VALUE'] = fom_dict[idx] / duration_dict[tech_params[0][3:6]]
        
    # Update VariableCosts with user-defined variable costs by storage technology
    df_var_cost = var_cost_base.copy()

    ''' Sets variable cost by taking the defined variable cost ($/MWh) divided by the storage 
    duration (=Storage Capacity (GWh)/Storage Power Rating (GW)) to mimic variable costs for storage
    as capacities are set in GWh/PJ terms.'''
    for idx, tech_params in tech_capacity_sto.items():
        df_var_cost.loc[df_var_cost['TECHNOLOGY'] == tech_params[0],
                   'VALUE'] = round(var_dict[idx] / duration_dict[tech_params[0][3:6]] / 3.6 , 4)
        
    return(df_max_cap_inv, 
           df_min_cap_inv, 
           df_res_cap,
           df_res_sto_cap,
           df_oar, 
           df_cap_cost, 
           df_cap_cost_sto,
           df_fix_cost,
           df_var_cost
           )