"""Function to calculate residual capacity for powerplant technologies."""

import pandas as pd
import itertools

from constants import (
    region_name,
    custom_nodes,
    years,
    new_iar_ccg,
    new_iar_ocg,
    new_iar_coa,
    new_iar_default,
)

from data import(
    createPwrTechs, 
    duplicatePlexosTechs,
    newIar
    )

from utils import apply_dtypes

def activity_master_start(df_gen_base, nodes_extra_list, duplicate_techs, mode_list):

    # ### Add input and output activity ratios

    # Create master table for activity ratios 
    if custom_nodes:
        node_list = list(df_gen_base['node_code'].unique()) + custom_nodes
    else:
        node_list = list(df_gen_base['node_code'].unique())
    # Add extra nodes which are not present in 2015 but will be by 2050
    for each_node in nodes_extra_list:
        if len(each_node) <= 6:
            node_list.append("".join(each_node.split('-')[1:]) + 'XX')
        else:
            node_list.append("".join(each_node.split('-')[1:]))

    master_fuel_list = list(df_gen_base['tech_code'].unique())
    master_fuel_list.append('CCS')

    df_ratios = pd.DataFrame(list(itertools.product(node_list,
                                                    master_fuel_list,
                                                    mode_list,
                                                    years)
                                  ),
                             columns = ['node_code', 'tech_code', 'MODE_OF_OPERATION', 'YEAR']
                             )

    df_ratios = createPwrTechs(df_ratios, duplicate_techs)

    return df_ratios

def activity_output_pwr(df_ratios, thermal_fuel_list_oar):

    # #### OutputActivityRatio - Power Generation Technologies
    df_oar = df_ratios.copy()
    mask = df_oar['TECHNOLOGY'].apply(lambda x: x[3:6] in thermal_fuel_list_oar)
    df_oar['FUEL'] = 0
    df_oar['FUEL'][mask] = 1
    df_oar = df_oar.loc[~((df_oar['MODE_OF_OPERATION'] > 1) &
                          (df_oar['FUEL'] == 0))]
    df_oar['FUEL'] = ('ELC' + 
                      df_oar['TECHNOLOGY'].str[6:11] + 
                      '01'
                     )
    df_oar['VALUE'] = 1
    
    # Add 'REGION' column and fill 'GLOBAL' throughout
    df_oar['REGION'] = region_name
    
    # Select columns for final output table
    df_oar_base = df_oar[['REGION', 
                     'TECHNOLOGY',
                     'FUEL',                  
                     'MODE_OF_OPERATION',
                     'YEAR', 
                     'VALUE',]]
    
    return df_oar, df_oar_base

def activity_input_pwr(df_oar, thermal_fuel_list_iar, renewables_list, 
                       df_eff_node, df_eff_tech):

    # #### InputActivityRatio - Power Generation Technologies
    # Copy OAR table with all columns to IAR
    df_iar = df_oar.copy()

    df_iar['FUEL'] = 0

    # Deal with GAS techs first...  OCG and CCG
    # OCG Mode 1: Domestic GAS
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(['OCG'])),
               'FUEL'] = 'GAS'+df_iar['TECHNOLOGY'].str[6:9]
    # OCG Mode 2: International GAS
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 2) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(['OCG'])),
               'FUEL'] = 'GASINT'

    # CCG Mode 1: Domestic GAS
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(['CCG'])),
               'FUEL'] = 'GAS'+df_iar['TECHNOLOGY'].str[6:9]

    # CCG Mode 2: International GAS
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 2) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(['CCG'])),
               'FUEL'] = 'GASINT'
    
    # CCS Mode 1: Domestic COA
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(['CCS'])),
               'FUEL'] = 'COA'+df_iar['TECHNOLOGY'].str[6:9]

    # CCS Mode 2: International COA
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 2) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(['CCS'])),
               'FUEL'] = 'COAINT'

    # For non-GAS thermal fuels, domestic fuel input by country in mode 1 and 
    # 'international' fuel input in mode 2
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(thermal_fuel_list_iar)),
               'FUEL'] = df_iar['TECHNOLOGY'].str[3:9]

    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 2) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(thermal_fuel_list_iar)),
               'FUEL'] = df_iar['TECHNOLOGY'].str[3:6] + 'INT'

    # For renewable fuels, input by node in mode 1
    df_iar.loc[(df_iar['MODE_OF_OPERATION'] == 1) &
               (df_iar['TECHNOLOGY'].str[3:6].isin(renewables_list)),
               'FUEL'] = df_iar['TECHNOLOGY'].str[3:11]

    # Remove mode 2 when not used
    df_iar = df_iar.loc[df_iar['FUEL'] != 0]

    # Join efficiency columns: one with node and technology average, and the
    # other with technology average
    df_iar = df_iar.join(df_eff_node.set_index(['tech_code', 'node_code']), 
                         on=['tech_code', 'node_code'])

    df_iar = df_iar.join(df_eff_tech.set_index('tech_code'), 
                         on='tech_code')

    # When available, choose node and technology average. Else, 
    # choose technology average
    df_iar['VALUE'] = df_iar['node_average_iar']
    
    df_iar.loc[df_iar['TECHNOLOGY'].str.startswith('PWRCCS'),
               'tech_average_iar'] = 3
    
    df_iar.loc[df_iar['VALUE'].isna(),
               'VALUE'] = df_iar['tech_average_iar']
    df_iar.drop_duplicates(inplace=True)

    # Add 'REGION' column and fill 'GLOBAL' throughout
    df_iar['REGION'] = region_name

    # Select columns for final output table
    df_iar_base = df_iar[['REGION', 
                           'TECHNOLOGY',
                           'FUEL',  
                           'MODE_OF_OPERATION',
                           'YEAR', 
                           'VALUE',]]
    
    return df_iar_base

def activity_upstream(df_iar_base, renewables_list, thermal_fuels_list_mining):

    # #### OutputActivityRatios - Upstream
    
    # We have to create a technology to produce every fuel that is input into any of the power technologies:
    
    df_oar_upstream = df_iar_base.copy()
    
    # All mining and resource technologies have an OAR of 1...
    df_oar_upstream['VALUE'] = 1
    
    # Renewables - set the technology as RNW + FUEL
    df_oar_upstream.loc[df_oar_upstream['FUEL'].str[0:3].isin(renewables_list),
               'TECHNOLOGY'] = 'RNW'+df_oar_upstream['FUEL']
    
    # If the fuel is a thermal fuel, we need to create the OAR for the mining technology... BUT NOT FOR THE INT FUELS...
    df_oar_upstream.loc[df_oar_upstream['FUEL'].str[0:3].isin(thermal_fuels_list_mining) & 
                        ~(df_oar_upstream['FUEL'].str[3:6] == "INT"),
               'TECHNOLOGY'] = 'MIN'+df_oar_upstream['FUEL']
    
    # Above should get all the outputs for the MIN technologies, but we need to adjust the mode 2 ones to just the fuel code (rather than MINCOAINT)
    df_oar_upstream.loc[df_oar_upstream['MODE_OF_OPERATION']==2,
               'TECHNOLOGY'] = 'MIN'+df_oar_upstream['FUEL'].str[0:3]+df_oar_upstream['TECHNOLOGY'].str[6:9]
    df_oar_upstream.loc[df_oar_upstream['MODE_OF_OPERATION']==2,
               'FUEL'] = df_oar_upstream['FUEL'].str[0:3]

    # Now remove the duplicate fuels that the above created (because there's now a COA for each country, not each region, and GAS is repeated twice for each region as well):
    df_oar_upstream.drop_duplicates(keep='first',inplace=True)
    
    # Now we have to create the MINXXXINT technologies.  They are all based on the MODE_OF_OPERATION == 2:
    df_oar_int = pd.DataFrame(df_oar_upstream.loc[df_oar_upstream['MODE_OF_OPERATION'] == 2, :])
    
    # At this point we should have only the internationally traded fuels since they're all mode 2.  So we can make the tech MINXXXINT and that's that.
    df_oar_int['TECHNOLOGY'] = 'MIN'+df_oar_int['FUEL']+'INT'
    # And rename the fuel to XXXINT
    df_oar_int['FUEL'] = df_oar_int['FUEL']+'INT'
    df_oar_int['MODE_OF_OPERATION'] = 1  # This is probably not strictly necessary as long as they're always the same in and out...
    
    # and de-duplicate this list:
    df_oar_int.drop_duplicates(keep='first',inplace=True)

    # #### Input Activity Ratios - Upstream
    
    # All we need to do is take in the thermal fuels for the MINXXXINT technologies.  This already exists as df_oar_int with the XXINT fuel so we can simply copy that:
    df_iar_int = df_oar_int.copy()
    df_iar_int['FUEL'] = df_iar_int['FUEL'].str[0:3]
    
    return df_oar_upstream, df_oar_int

def activity_master_end(df_oar_base, df_oar_upstream, df_oar_int, df_oar_trn, 
                        df_int_trn_oar, df_iar_final, df_iar_trn, df_int_trn_iar,
                        duplicate_techs):

    # #### Output IAR and OAR

    # Combine the pieces from above and output to csv:
    
    df_oar_final = pd.concat(
        [
            df_oar_base, # If we want to split transmission to different script might have to not used df_oar_final as input to this function.
            df_oar_upstream,
            df_oar_int,
            df_oar_trn,
            df_int_trn_oar,
        ]
    ).dropna()
    
    # Select columns for final output table
    df_oar_final = df_oar_final[['REGION', 
                                 'TECHNOLOGY',
                                 'FUEL',  
                                 'MODE_OF_OPERATION',
                                 'YEAR', 
                                 'VALUE',]]
    
    df_iar_final = pd.concat(
        [
            df_iar_final,
            df_iar_trn,
            df_int_trn_iar,
        ]
    ).dropna()
    
    # Select columns for final output table
    df_iar_final = df_iar_final[['REGION', 
                                 'TECHNOLOGY',
                                 'FUEL',  
                                 'MODE_OF_OPERATION',
                                 'YEAR', 
                                 'VALUE']]
    
    # Add iar for techs not using PLEXOS values 
    df_iar_newTechs = duplicatePlexosTechs(df_iar_final, duplicate_techs)
    for duplicate_tech in duplicate_techs:
        df_new_iar = newIar(df_iar_newTechs, duplicate_tech,new_iar_ccg, 
                   new_iar_ocg, new_iar_coa, new_iar_default)
        df_iar_final = pd.concat([df_iar_final, df_new_iar])
    
    # Add oar for techs not using PLEXOS values 
    df_oar_newTechs = duplicatePlexosTechs(df_oar_final, duplicate_techs)
    df_oar_final = pd.concat([df_oar_final, df_oar_newTechs])

    df_oar_final.drop_duplicates(subset=['REGION', 
                                         'TECHNOLOGY',
                                         'FUEL',  
                                         'MODE_OF_OPERATION',
                                         'YEAR'],
                                 keep='last',
                                 inplace=True)
    
    df_iar_final = apply_dtypes(df_iar_final, "InputActivityRatio")
    
    return df_oar_final, df_iar_final

def capact(df_oar_final):

    # Create CapacityToActivityUnit csv
    df_capact_final = df_oar_final[['REGION',
                                    'TECHNOLOGY'
                                    ]]
    df_capact_final = df_capact_final.drop_duplicates()
    df_capact_final = (df_capact_final
                       .loc[(df_capact_final['TECHNOLOGY']
                             .str.startswith('PWR')
                             ) | 
                            (df_capact_final['TECHNOLOGY']
                             .str.contains('TRN')
                             )
                            ]
                       )

    df_capact_final['VALUE'] = 31.536
    df_capact_final.drop_duplicates(inplace=True)

    return df_capact_final