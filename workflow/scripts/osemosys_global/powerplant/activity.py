"""Function to calculate activity for powerplant technologies."""

import pandas as pd
import itertools

from constants import (
    NODES_EXTRA_LIST,
    NEW_IAR_CCG,
    NEW_IAR_OCG,
    NEW_IAR_COA,
    NEW_IAR_DEFAULT,
    THERMAL_FUEL_LIST_OAR,
    THERMAL_FUEL_LIST_IAR,
    THERMAL_FUEL_LIST_MINING,
)

from data import(
    create_pwr_techs, 
    duplicate_plexos_techs,
    new_iar,
    get_years
    )

from utils import apply_dtypes

def activity_master_start(df_gen_base, duplicate_techs, mode_list,
                          custom_nodes, start_year, end_year):

    # ### Add input and output activity ratios

    # Create master table for activity ratios 
    if custom_nodes:
        node_list = list(df_gen_base['node_code'].unique()) + custom_nodes
    else:
        node_list = list(df_gen_base['node_code'].unique())
    # Add extra nodes which are not present in 2015 but will be by 2050
    for each_node in NODES_EXTRA_LIST:
        if len(each_node) <= 6:
            node_list.append("".join(each_node.split('-')[1:]) + 'XX')
        else:
            node_list.append("".join(each_node.split('-')[1:]))

    master_fuel_list = list(df_gen_base['tech_code'].unique())
    master_fuel_list.append('CCS')

    df_ratios = pd.DataFrame(list(itertools.product(node_list,
                                                    master_fuel_list,
                                                    mode_list,
                                                    get_years(start_year, end_year))
                                  ),
                             columns = ['node_code', 'tech_code', 'MODE_OF_OPERATION', 'YEAR']
                             )

    df_ratios = create_pwr_techs(df_ratios, duplicate_techs)

    return df_ratios

def activity_output_pwr(df_ratios, region_name):

    # #### OutputActivityRatio - Power Generation Technologies
    df_pwr_oar_base = df_ratios.copy()
    mask = df_pwr_oar_base['TECHNOLOGY'].apply(lambda x: x[3:6] in THERMAL_FUEL_LIST_OAR)
    df_pwr_oar_base['FUEL'] = 0
    df_pwr_oar_base.loc[mask, "FUEL"] = 1
    
    df_pwr_oar_base = df_pwr_oar_base.loc[~((df_pwr_oar_base['MODE_OF_OPERATION'] > 1) &
                          (df_pwr_oar_base['FUEL'] == 0))]
    df_pwr_oar_base['FUEL'] = ('ELC' + 
                      df_pwr_oar_base['TECHNOLOGY'].str[6:11] + 
                      '01'
                     )
    df_pwr_oar_base['VALUE'] = 1
    
    # Add 'REGION' column and fill 'GLOBAL' throughout
    df_pwr_oar_base['REGION'] = region_name
    
    # Select columns for final output table
    df_pwr_oar_final = df_pwr_oar_base[['REGION', 
                     'TECHNOLOGY',
                     'FUEL',                  
                     'MODE_OF_OPERATION',
                     'YEAR', 
                     'VALUE',]]
    
    # Copy OAR table with all columns to IAR
    df_pwr_iar_base = df_pwr_oar_base.copy()
    
    return df_pwr_oar_final, df_pwr_iar_base

def activity_input_pwr(df_pwr_iar_base, renewables_list, df_eff_node, 
                       df_eff_tech, region_name):

    # #### InputActivityRatio - Power Generation Technologies
    df_pwr_iar_base['FUEL'] = 0
    df_pwr_iar_base['FUEL'] = df_pwr_iar_base['FUEL'].astype(str)

    # Deal with GAS techs first...  OCG and CCG
    # OCG Mode 1: Domestic GAS
    df_pwr_iar_base.loc[(df_pwr_iar_base['MODE_OF_OPERATION'] == 1) &
               (df_pwr_iar_base['TECHNOLOGY'].str[3:6].isin(['OCG'])),
               'FUEL'] = 'GAS'+df_pwr_iar_base['TECHNOLOGY'].str[6:9]
    # OCG Mode 2: International GAS
    df_pwr_iar_base.loc[(df_pwr_iar_base['MODE_OF_OPERATION'] == 2) &
               (df_pwr_iar_base['TECHNOLOGY'].str[3:6].isin(['OCG'])),
               'FUEL'] = 'GASINT'

    # CCG Mode 1: Domestic GAS
    df_pwr_iar_base.loc[(df_pwr_iar_base['MODE_OF_OPERATION'] == 1) &
               (df_pwr_iar_base['TECHNOLOGY'].str[3:6].isin(['CCG'])),
               'FUEL'] = 'GAS'+df_pwr_iar_base['TECHNOLOGY'].str[6:9]

    # CCG Mode 2: International GAS
    df_pwr_iar_base.loc[(df_pwr_iar_base['MODE_OF_OPERATION'] == 2) &
               (df_pwr_iar_base['TECHNOLOGY'].str[3:6].isin(['CCG'])),
               'FUEL'] = 'GASINT'
    
    # CCS Mode 1: Domestic COA
    df_pwr_iar_base.loc[(df_pwr_iar_base['MODE_OF_OPERATION'] == 1) &
               (df_pwr_iar_base['TECHNOLOGY'].str[3:6].isin(['CCS'])),
               'FUEL'] = 'COA'+df_pwr_iar_base['TECHNOLOGY'].str[6:9]

    # CCS Mode 2: International COA
    df_pwr_iar_base.loc[(df_pwr_iar_base['MODE_OF_OPERATION'] == 2) &
               (df_pwr_iar_base['TECHNOLOGY'].str[3:6].isin(['CCS'])),
               'FUEL'] = 'COAINT'

    # For non-GAS thermal fuels, domestic fuel input by country in mode 1 and 
    # 'international' fuel input in mode 2
    df_pwr_iar_base.loc[(df_pwr_iar_base['MODE_OF_OPERATION'] == 1) &
               (df_pwr_iar_base['TECHNOLOGY'].str[3:6].isin(THERMAL_FUEL_LIST_IAR)),
               'FUEL'] = df_pwr_iar_base['TECHNOLOGY'].str[3:9]

    df_pwr_iar_base.loc[(df_pwr_iar_base['MODE_OF_OPERATION'] == 2) &
               (df_pwr_iar_base['TECHNOLOGY'].str[3:6].isin(THERMAL_FUEL_LIST_IAR)),
               'FUEL'] = df_pwr_iar_base['TECHNOLOGY'].str[3:6] + 'INT'

    # For renewable fuels, input by node in mode 1
    df_pwr_iar_base.loc[(df_pwr_iar_base['MODE_OF_OPERATION'] == 1) &
               (df_pwr_iar_base['TECHNOLOGY'].str[3:6].isin(renewables_list)),
               'FUEL'] = df_pwr_iar_base['TECHNOLOGY'].str[3:11]

    # Remove mode 2 when not used
    df_pwr_iar_base = df_pwr_iar_base.loc[df_pwr_iar_base['FUEL'] != 0]

    # Join efficiency columns: one with node and technology average, and the
    # other with technology average
    df_pwr_iar_base = df_pwr_iar_base.join(df_eff_node.set_index(['tech_code', 'node_code']), 
                         on=['tech_code', 'node_code'])

    df_pwr_iar_base = df_pwr_iar_base.join(df_eff_tech.set_index('tech_code'), 
                         on='tech_code')

    # When available, choose node and technology average. Else, 
    # choose technology average
    df_pwr_iar_base['VALUE'] = df_pwr_iar_base['node_average_iar']
    
    df_pwr_iar_base.loc[df_pwr_iar_base['TECHNOLOGY'].str.startswith('PWRCCS'),
               'tech_average_iar'] = 3
    
    df_pwr_iar_base.loc[df_pwr_iar_base['VALUE'].isna(),
               'VALUE'] = df_pwr_iar_base['tech_average_iar']
    df_pwr_iar_base.drop_duplicates(inplace=True)

    # Add 'REGION' column and fill 'GLOBAL' throughout
    df_pwr_iar_base['REGION'] = region_name

    # Select columns for final output table
    df_pwr_iar_final = df_pwr_iar_base[['REGION', 
                                        'TECHNOLOGY',
                                        'FUEL',  
                                        'MODE_OF_OPERATION',
                                        'YEAR', 
                                        'VALUE',]]
    
    return df_pwr_iar_final

def activity_upstream(df_pwr_iar_final, renewables_list):

    # #### OutputActivityRatios - Upstream
    
    # We have to create a technology to produce every fuel that is input into any of the power technologies:
    
    df_oar_upstream = df_pwr_iar_final.copy()
    
    # All mining and resource technologies have an OAR of 1...
    df_oar_upstream['VALUE'] = 1
    
    # Renewables - set the technology as RNW + FUEL
    df_oar_upstream.loc[df_oar_upstream['FUEL'].str[0:3].isin(renewables_list),
               'TECHNOLOGY'] = 'RNW'+df_oar_upstream['FUEL']
    
    # If the fuel is a thermal fuel, we need to create the OAR for the mining technology... BUT NOT FOR THE INT FUELS...
    df_oar_upstream.loc[df_oar_upstream['FUEL'].str[0:3].isin(THERMAL_FUEL_LIST_MINING) & 
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

def activity_master_end(df_pwr_oar_final, df_oar_upstream, df_oar_int, 
                        df_pwr_iar_final, duplicate_techs):

    # #### Output IAR and OAR

    # Combine the pieces from above and output to csv:
    
    df_oar_final = pd.concat(
        [
            df_pwr_oar_final,
            df_oar_upstream,
            df_oar_int,
        ]
    ).dropna()
    
    # Select columns for final output table
    df_oar_final = df_oar_final[['REGION', 
                                 'TECHNOLOGY',
                                 'FUEL',  
                                 'MODE_OF_OPERATION',
                                 'YEAR', 
                                 'VALUE',]]
    
    # Select columns for final output table
    df_iar_final = df_pwr_iar_final[['REGION', 
                                 'TECHNOLOGY',
                                 'FUEL',  
                                 'MODE_OF_OPERATION',
                                 'YEAR', 
                                 'VALUE']]
    
    # Add iar for techs not using PLEXOS values 
    df_iar_newTechs = duplicate_plexos_techs(df_iar_final, duplicate_techs)
    for duplicate_tech in duplicate_techs:
        df_new_iar = new_iar(df_iar_newTechs, duplicate_tech,NEW_IAR_CCG, 
                   NEW_IAR_OCG, NEW_IAR_COA, NEW_IAR_DEFAULT)
        df_iar_final = pd.concat([df_iar_final, df_new_iar])
    
    # Add oar for techs not using PLEXOS values 
    df_oar_newTechs = duplicate_plexos_techs(df_oar_final, duplicate_techs)
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
    df_capact_final = df_capact_final.loc[df_capact_final['TECHNOLOGY'
                                                          ].str.startswith('PWR')]

    df_capact_final['VALUE'] = 31.536
    df_capact_final.drop_duplicates(inplace=True)

    return df_capact_final

def set_annual_activity_upper_limit(fuel_limit, start_year, 
                                    end_year, region_name):
    
    years = get_years(start_year, end_year)
    
    mf_df = fuel_limit.copy()
    mf_df["TECHNOLOGY"] = "MIN" + mf_df["FUEL"] + mf_df["COUNTRY"]
    mf_df = mf_df[["TECHNOLOGY", "YEAR", "VALUE"]]

    tech_list = mf_df["TECHNOLOGY"].unique()
    mf_df_final = pd.DataFrame(
        list(itertools.product(tech_list, years)), columns=["TECHNOLOGY", "YEAR"]
    )
    mf_df_final = pd.merge(mf_df_final, mf_df, how="left", on=["TECHNOLOGY", "YEAR"])
    mf_df_final["VALUE"] = mf_df_final["VALUE"].astype(float)
    for each_tech in tech_list:
        mf_df_final.loc[mf_df_final["TECHNOLOGY"].isin([each_tech]), "VALUE"] = (
            mf_df_final.loc[mf_df_final["TECHNOLOGY"].isin([each_tech]), "VALUE"]
            .interpolate()
            .round(0)
        )

    mf_df_final["REGION"] = region_name
    mf_df_final = mf_df_final[["REGION", "TECHNOLOGY", "YEAR", "VALUE"]]
    mf_df_final.dropna(inplace=True)

    return mf_df_final

def set_model_period_activity_upper_limit(tech_set, region_name):

    # Model Period Activity Upper Limit for 'MINCOA***01'
    min_tech_df = tech_set.copy()
    min_tech = [
        x
        for x in min_tech_df["VALUE"].unique()
        if x.startswith("MINCOA")
        if x.endswith("01")
    ]
    min_tech_df_final = pd.DataFrame(columns=["REGION", "TECHNOLOGY", "VALUE"])
    min_tech_df_final["TECHNOLOGY"] = min_tech
    min_tech_df_final["REGION"] = region_name
    min_tech_df_final["VALUE"] = 0

    return min_tech_df_final