"""Function to integrate user defined capacities."""

import pandas as pd
import itertools

from constants import (
    region_name,
    years
)

def set_user_defined_capacity(tech_capacity, op_life_dict, df_tech_set, 
                      df_min_cap_invest, df_max_cap_invest, df_res_cap,
                      df_iar_final, df_oar_final, fuel_set, op_life_base, cap_act_base,
                      cap_cost_base, df_iar_custom_val, df_oar_custom_val):
    """User-defined capacities are used when a specific technology must be 
    invested, for a given year and capacity. This is applied through the 
    parameter 'TotalAnnualMinCapacityInvestment'. 

    Args:
        region: From config file (e.g. 'GLOBAL')
        tech_capacity: User-defined capacity in config file 
                       (e.g. TRNAGOXXCODXX: [5, 2030])

    Returns
        None
    """
    techCapacity = []
    tech_capacity_dict = {}
    first_year_dict = {}
    build_rate_dict = {}
    capex_dict = {}
    build_year_dict = {}
    
    if not tech_capacity is None:

        for tech, tech_params in tech_capacity.items():
            techCapacity.append([tech, tech_params[0], tech_params[1]])
            tech_capacity_dict[tech] = tech_params[2]
            build_year_dict[tech] = tech_params[1]
            first_year_dict[tech] = tech_params[3]
            build_rate_dict[tech] = tech_params[4]
            capex_dict[tech] = tech_params[5] # 
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
        
        df = pd.DataFrame(list(itertools.product(list(tech_capacity_df['TECHNOLOGY'].unique()),
                                                 years)
                               ),
                          columns = ['TECHNOLOGY', 
                                     'YEAR']
                          )
        df['REGION'] = region_name
        df = pd.merge(df, df_min_cap_inv,
                      how='left',
                      on=['REGION', 'TECHNOLOGY', 'YEAR'])
        df['FIRST_YEAR'] = df['TECHNOLOGY'].map(first_year_dict)
        df['BUILD_YEAR'] = df['TECHNOLOGY'].map(build_year_dict)
        df['MAX_BUILD'] = df['TECHNOLOGY'].map(build_rate_dict)

        # Fill VALUE with MAX_BUILD for YEAR >= FIRST_YEAR
        df.loc[(df['YEAR']>=df['FIRST_YEAR']) &
               (df['YEAR']>df['BUILD_YEAR']),
               'VALUE'] = df['MAX_BUILD']
        df.infer_objects().fillna(0,
                  inplace=True)
        max_cap_techs_df = df[['REGION',
                               'TECHNOLOGY',
                               'YEAR',
                               'VALUE']]

        #Replace VALUE with CAPEX
        df['VALUE'] = df['TECHNOLOGY'].map(capex_dict)
        
        # Append existing TotalAnnualMaxCapacityInvestment data with MAX_BUILD for TRN
        df_max_cap_inv = pd.concat([df_max_cap_invest, max_cap_techs_df]).drop_duplicates()

        # Print TotalAnnualMaxCapacityInvestment.csv with MAX_BUILD for TRN
        df_max_cap_inv.drop_duplicates(subset=['REGION', 
                                               'TECHNOLOGY',
                                               'YEAR'],
                                       keep='last',
                                       inplace=True)
        
        # For technologies with start year before model start year, add to ResidualCapacity
        df_res_cap_ud = df_min_cap_inv.copy().loc[df_min_cap_inv.copy()['YEAR'] < min(years)]
        df_res_cap_ud.rename(columns={'YEAR':'START_YEAR'},
                             inplace=True)
        df_res_cap_ud_final = pd.DataFrame(list(itertools.product(df_res_cap_ud['TECHNOLOGY'].unique(),
                                                                  years)
                                                ),
                                           columns = ['TECHNOLOGY',
                                                      'YEAR']
                                           )
        df_res_cap_ud_final = pd.merge(df_res_cap_ud_final,
                                       df_res_cap_ud,
                                       how='left',
                                       on=['TECHNOLOGY'])
        df_res_cap_ud_final['TECH'] = df_res_cap_ud_final['TECHNOLOGY'].str[3:6]
        df_res_cap_ud_final.loc[df_res_cap_ud_final['TECHNOLOGY'].str.contains('TRN'),
                                'TECH'] = 'TRN'
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

        
        df_res_cap = pd.concat([df_res_cap, df_res_cap_ud_final 
                                if not df_res_cap_ud_final.empty else None])
                
        # For technologies with start year at or after model start year, add to 
        # TotalAnnualMinCapacityInvestment      
        df_min_cap_inv = df_min_cap_inv.loc[df_min_cap_inv['YEAR'] >= min(years)]
        df_min_cap_inv.drop_duplicates(subset=['REGION', 
                                               'TECHNOLOGY',
                                               'YEAR'],
                                       keep='last',
                                       inplace=True)
        
        # Add IAR and OAR for custom technologies
        tech_list = list(tech_capacity_df['TECHNOLOGY'].unique())
        df_iar_custom = pd.DataFrame(list(itertools.product(tech_list,
                                                        [1, 2],
                                                        years)
                                      ),
                                 columns = ['TECHNOLOGY',
                                            'MODE_OF_OPERATION',
                                            'YEAR']
                                 )
        df_oar_custom = pd.DataFrame(list(itertools.product(tech_list,
                                                        [1, 2],
                                                        years)
                                      ),
                                 columns = ['TECHNOLOGY',
                                            'MODE_OF_OPERATION',
                                            'YEAR']
                                 )  
        # IAR in modes 1 and 2 are primary electricity commodity ('ELC*01') in 
        # node_from and node_to, respectively. 
        # OAR is the inverse of the above
        df_iar_custom.loc[df_iar_custom['MODE_OF_OPERATION']==1,
                                        'FUEL'] = ('ELC' + 
                                                   df_iar_custom['TECHNOLOGY'].str[3:8] + 
                                                   '01')
        df_iar_custom.loc[df_iar_custom['MODE_OF_OPERATION']==2,
                                        'FUEL'] = ('ELC' + 
                                                   df_iar_custom['TECHNOLOGY'].str[8:13] + 
                                                   '01')
        df_oar_custom.loc[df_iar_custom['MODE_OF_OPERATION']==1,
                                        'FUEL'] = ('ELC' + 
                                                   df_oar_custom['TECHNOLOGY'].str[8:13] + 
                                                   '01')
        df_oar_custom.loc[df_iar_custom['MODE_OF_OPERATION']==2,
                                        'FUEL'] = ('ELC' + 
                                                   df_oar_custom['TECHNOLOGY'].str[3:8] + 
                                                   '01')
        df_iar_custom['VALUE'] = df_iar_custom_val
        df_oar_custom['VALUE'] = df_oar_custom_val
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
        
        df_iar = pd.concat([df_iar_final, df_iar_custom])
        df_oar = pd.concat([df_oar_final, df_oar_custom])
        
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

        op_life_custom.loc[op_life_custom['TECHNOLOGY'].str.contains('TRN'),
                           'VALUE'] = op_life_dict.get('TRN')
        op_life_custom['REGION'] = region_name
        op_life_custom = op_life_custom[['REGION',
                                         'TECHNOLOGY',
                                         'VALUE']]

        op_life = pd.concat([op_life_base, op_life_custom])
        op_life.drop_duplicates(subset=['REGION', 
                                        'TECHNOLOGY'],
                                keep='last',
                                inplace=True)

        cap_act_custom = pd.DataFrame({'TECHNOLOGY': tech_list})
        cap_act_custom.loc[cap_act_custom['TECHNOLOGY'].str.contains('TRN'),
                           'VALUE'] = 31.536
        cap_act_custom.loc[cap_act_custom['TECHNOLOGY'].str.contains('PWR'),
                           'VALUE'] = 31.536
        cap_act_custom['REGION'] = region_name
        cap_act_custom = cap_act_custom[['REGION',
                                         'TECHNOLOGY',
                                         'VALUE']]
        cap_act = pd.concat([cap_act_base, cap_act_custom])
        cap_act.drop_duplicates(inplace=True)
        
        tech_list = list(tech_capacity_df['TECHNOLOGY'].unique())
        cap_cost_trn = pd.DataFrame(list(itertools.product(tech_list,
                                                           years)),
                                    columns = ['TECHNOLOGY',
                                               'YEAR'])

        # Update CapitalCost with user-defined costs by transmission line
        for each_trn in tech_list:
            cap_cost_trn.loc[cap_cost_trn['TECHNOLOGY'].str.startswith(each_trn),
                             'VALUE'] = capex_dict[each_trn]

        cap_cost_trn.loc[cap_cost_trn['TECHNOLOGY'].str.contains('PWRTRN'),
                         'VALUE'] = 300
        cap_cost_trn['REGION'] = region_name
        cap_cost_trn = cap_cost_trn[['REGION',
                                     'TECHNOLOGY',
                                     'YEAR',
                                     'VALUE']]
        cap_cost = pd.concat([cap_cost_base, cap_cost_trn])
        cap_cost.drop_duplicates(subset=['REGION', 'TECHNOLOGY', 'YEAR'],
                                 keep="last",
                                 inplace=True)
        
        return(df_tech_set, fuel_set, df_max_cap_inv, df_min_cap_inv, 
               df_res_cap, df_iar, df_oar, op_life, cap_act, cap_cost)
        