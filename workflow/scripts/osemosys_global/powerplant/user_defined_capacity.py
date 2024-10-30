"""Function to integrate user defined capacities."""

import pandas as pd
import itertools

from data import get_years

def set_user_defined_capacity(tech_capacity, op_life_dict, df_tech_set, 
                      df_min_cap_invest, df_max_cap_invest, df_res_cap,
                      op_life_base, cap_act_base, cap_cost_base, df_iar_final,
                      df_oar_final, fuel_set, start_year, end_year, region_name,
                      renewables_list
                      ):
    
    techCapacity = []
    first_year_expansion_dict = {}
    build_rate_dict = {}
    capex_dict = {}
    build_year_dict = {}
    efficiency_dict = {}

    for tech, tech_params in tech_capacity.items():
        techCapacity.append([tech, tech_params[0], tech_params[1]])
        build_year_dict[tech] = tech_params[1]
        first_year_expansion_dict[tech] = tech_params[2]
        build_rate_dict[tech] = tech_params[3]
        capex_dict[tech] = tech_params[4]
        efficiency_dict[tech] = tech_params[5]       
    tech_capacity_df = pd.DataFrame(techCapacity,
                                    columns=['TECHNOLOGY', 'VALUE', 'YEAR'])
    tech_capacity_df['REGION'] = region_name
    tech_capacity_df = tech_capacity_df[['REGION', 'TECHNOLOGY', 'YEAR', 'VALUE']]

    for each_tech in list(tech_capacity_df['TECHNOLOGY'].unique()):
        if each_tech not in list(df_tech_set['VALUE']):
            df_tech_set = pd.concat([df_tech_set, pd.DataFrame({'VALUE':[each_tech]})])
            
    df_tech_set.drop_duplicates(inplace=True)

    df_min_cap_inv = pd.concat([df_min_cap_invest, tech_capacity_df])
    df_min_cap_inv.drop_duplicates(inplace=True)
    
    max_cap_techs_df = pd.DataFrame(list(itertools.product(list(tech_capacity_df['TECHNOLOGY'].unique()),
                                             get_years(start_year, end_year))
                           ),
                      columns = ['TECHNOLOGY', 
                                 'YEAR']
                      )
    max_cap_techs_df['REGION'] = region_name
    max_cap_techs_df = pd.merge(max_cap_techs_df, df_min_cap_inv,
                  how='left',
                  on=['REGION', 'TECHNOLOGY', 'YEAR'])
    max_cap_techs_df['FIRST_YEAR'] = max_cap_techs_df['TECHNOLOGY'].map(first_year_expansion_dict)
    max_cap_techs_df['BUILD_YEAR'] = max_cap_techs_df['TECHNOLOGY'].map(build_year_dict)
    max_cap_techs_df['MAX_BUILD'] = max_cap_techs_df['TECHNOLOGY'].map(build_rate_dict)

    # Fill VALUE with MAX_BUILD for YEAR >= FIRST_YEAR
    max_cap_techs_df.loc[(max_cap_techs_df['YEAR']>=max_cap_techs_df['FIRST_YEAR']) &
           (max_cap_techs_df['YEAR']>max_cap_techs_df['BUILD_YEAR']),
           'VALUE'] = max_cap_techs_df['MAX_BUILD']
    max_cap_techs_df['VALUE'] = max_cap_techs_df['VALUE'].infer_objects(copy=False).fillna(0)
    
    max_cap_techs_df = max_cap_techs_df[['REGION',
                                         'TECHNOLOGY',
                                         'YEAR',
                                         'VALUE']]
    
    # Append existing TotalAnnualMaxCapacityInvestment data with MAX_BUILD
    df_max_cap_inv = pd.concat([df_max_cap_invest, max_cap_techs_df]).drop_duplicates()

    # Print TotalAnnualMaxCapacityInvestment.csv with MAX_BUILD
    df_max_cap_inv.drop_duplicates(subset=['REGION', 
                                           'TECHNOLOGY',
                                           'YEAR'],
                                   keep='last',
                                   inplace=True)

    # For technologies with start year before model start year, add to ResidualCapacity
    df_res_cap_ud = df_min_cap_inv.copy().loc[df_min_cap_inv.copy()['YEAR'] < min(get_years(start_year, end_year))]
    df_res_cap_ud.rename(columns={'YEAR':'start_year'},
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
    df_res_cap_ud_final['end_year'] = (df_res_cap_ud_final['OP_LIFE'] 
                                       + df_res_cap_ud_final['start_year'])
    df_res_cap_ud_final = df_res_cap_ud_final.loc[df_res_cap_ud_final['YEAR'] 
                                                  >= df_res_cap_ud_final['start_year']]
    df_res_cap_ud_final = df_res_cap_ud_final.loc[df_res_cap_ud_final['YEAR'] 
                                                  <= df_res_cap_ud_final['end_year']]
    df_res_cap_ud_final['REGION'] = region_name
    df_res_cap_ud_final = df_res_cap_ud_final[['REGION',
                                               'TECHNOLOGY',
                                               'YEAR',
                                               'VALUE']]

    if not df_res_cap_ud_final.empty:
        df_res_cap = pd.concat([df_res_cap, df_res_cap_ud_final])
        df_res_cap = df_res_cap.groupby(by=["REGION", "TECHNOLOGY", "YEAR"], as_index=False).sum()
            
    # For technologies with start year at or after model start year, add to 
    # TotalAnnualMinCapacityInvestment      
    df_min_cap_inv = df_min_cap_inv.loc[df_min_cap_inv['YEAR'] >= min(get_years(start_year, end_year))]
    df_min_cap_inv.drop_duplicates(subset=['REGION', 
                                           'TECHNOLOGY',
                                           'YEAR'],
                                   keep='last',
                                   inplace=True)
    
    tech_list = list(tech_capacity_df['TECHNOLOGY'].unique())
   
    df_iar_custom = pd.DataFrame(list(itertools.product(tech_list,
                                                    [1, 2],
                                                    get_years(start_year, end_year))
                                  ),
                             columns = ['TECHNOLOGY',
                                        'MODE_OF_OPERATION',
                                        'YEAR']
                             )
    df_oar_custom = pd.DataFrame(list(itertools.product(tech_list,
                                                    [1, 2],
                                                    get_years(start_year, end_year))
                                  ),
                             columns = ['TECHNOLOGY',
                                        'MODE_OF_OPERATION',
                                        'YEAR']
                             )  

    for each_tech in tech_list:
        if each_tech[3:6] in renewables_list:

            df_iar_custom.loc[(df_iar_custom['TECHNOLOGY'] == each_tech) & 
                              (df_iar_custom['MODE_OF_OPERATION']==1),'FUEL'] = (
                                  df_iar_custom['TECHNOLOGY'].str[3:11])           
            
            df_oar_custom.loc[(df_oar_custom['TECHNOLOGY'] == each_tech) & 
                              (df_oar_custom['MODE_OF_OPERATION']==1),'FUEL'] = (
                'ELC' + df_oar_custom['TECHNOLOGY'].str[6:11] + '01')
            
        else:

            df_iar_custom.loc[(df_iar_custom['TECHNOLOGY'] == each_tech) & 
                              (df_iar_custom['MODE_OF_OPERATION']==1),'FUEL'] = (
                                  df_iar_custom['TECHNOLOGY'].str[3:9])
                                  
            df_iar_custom.loc[(df_iar_custom['TECHNOLOGY'] == each_tech) & 
                              (df_iar_custom['MODE_OF_OPERATION']==2),'FUEL'] = (
                                  df_iar_custom['TECHNOLOGY'].str[3:6] + "INT")
            
            df_oar_custom.loc[(df_oar_custom['TECHNOLOGY'] == each_tech) & 
                              (df_oar_custom['MODE_OF_OPERATION']==1),'FUEL'] = (
                'ELC' + df_oar_custom['TECHNOLOGY'].str[6:11] + '01')
            df_oar_custom.loc[(df_oar_custom['TECHNOLOGY'] == each_tech) & 
                              (df_oar_custom['MODE_OF_OPERATION']==2),'FUEL'] = (
                'ELC' + df_oar_custom['TECHNOLOGY'].str[6:11] + '01')
            
    for each_tech in tech_list:
        df_iar_custom.loc[df_iar_custom['TECHNOLOGY'] == each_tech,
                         'VALUE'] = round(1 / (efficiency_dict[each_tech] / 100), 3)

    df_oar_custom['VALUE'] = 1
    df_iar_custom['REGION'] = region_name
    df_oar_custom['REGION'] = region_name

    df_iar_custom = df_iar_custom[['REGION', 
                                   'TECHNOLOGY',
                                   'FUEL', 
                                   'MODE_OF_OPERATION',
                                   'YEAR', 
                                   'VALUE',]]
    df_oar_custom = df_oar_custom[['REGION', 
                                   'TECHNOLOGY',
                                   'FUEL', 
                                   'MODE_OF_OPERATION',
                                   'YEAR', 
                                   'VALUE',]]
    
    # Drop mode 2 for renewable techs
    df_iar_custom = df_iar_custom.loc[df_iar_custom['FUEL'] != 0]
    df_oar_custom = df_oar_custom.loc[df_oar_custom['FUEL'] != 0]
    
    df_iar = pd.concat([df_iar_final, df_iar_custom]).dropna()
    df_oar = pd.concat([df_oar_final, df_oar_custom]).dropna()
    
    df_iar.drop_duplicates(subset=['REGION', 
                                     'TECHNOLOGY',
                                     'FUEL',  
                                     'MODE_OF_OPERATION',
                                     'YEAR'],
                           keep='last',
                           inplace=True)

    df_oar.drop_duplicates(subset=['REGION', 
                                     'TECHNOLOGY',
                                     'FUEL',  
                                     'MODE_OF_OPERATION',
                                     'YEAR'],
                             keep='last',
                             inplace=True)

    # Add new fuels to FUEL set, if not already present
    fuel_list = []
    fuel_list = list(df_iar_custom['FUEL'].unique()) + list(df_oar_custom['FUEL'].unique())
    fuel_list = list(set(fuel_list))
    
    for each_fuel in fuel_list:
        if each_fuel not in list(fuel_set['VALUE']):
            fuel_set = pd.concat([fuel_set, pd.DataFrame({'VALUE':[each_fuel]})])          
    
    op_life_custom = pd.DataFrame({'TECHNOLOGY': tech_list})
    
    for each_tech in op_life_custom['TECHNOLOGY']:
        op_life_custom.loc[op_life_custom['TECHNOLOGY'] == each_tech,
                           'VALUE'] = op_life_dict.get(each_tech[3:6])

    op_life_custom['REGION'] = region_name
    op_life_custom = op_life_custom[['REGION',
                                     'TECHNOLOGY',
                                     'VALUE']]

    op_life = pd.concat([op_life_base, op_life_custom])
    op_life.drop_duplicates(subset=['REGION', 
                                    'TECHNOLOGY'],
                            keep='last',
                            inplace=True)

    # Add CapacityToActivityUnit for custom technologies
    cap_act_custom = pd.DataFrame({'TECHNOLOGY': tech_list})
    cap_act_custom.loc[cap_act_custom['TECHNOLOGY'].str.contains('PWR'),
                       'VALUE'] = 31.536
    cap_act_custom['REGION'] = region_name
    cap_act_custom = cap_act_custom[['REGION',
                                     'TECHNOLOGY',
                                     'VALUE']]
    cap_act = pd.concat([cap_act_base, cap_act_custom])
    cap_act.drop_duplicates(inplace=True)
    
    # Update CapitalCost with user-defined costs
    tech_list = list(tech_capacity_df['TECHNOLOGY'].unique())
    cap_cost = pd.DataFrame(list(itertools.product(tech_list,
                                                       get_years(start_year, end_year))),
                                columns = ['TECHNOLOGY',
                                           'YEAR'])

    for each_tech in tech_list:
        cap_cost.loc[cap_cost['TECHNOLOGY'].str.startswith(each_tech),
                         'VALUE'] = capex_dict[each_tech]

    cap_cost['REGION'] = region_name
    cap_cost = cap_cost[['REGION',
                                 'TECHNOLOGY',
                                 'YEAR',
                                 'VALUE']]
    cap_cost = pd.concat([cap_cost_base, cap_cost])
    cap_cost.drop_duplicates(subset=['REGION', 'TECHNOLOGY', 'YEAR'],
                             keep="last",
                             inplace=True)
    
    return(df_tech_set, df_max_cap_inv, df_min_cap_inv, 
           df_res_cap, df_iar, df_oar, fuel_set, op_life, 
           cap_act, cap_cost)