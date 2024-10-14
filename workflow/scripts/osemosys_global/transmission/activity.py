"""Functions to calaculate activity related to transmission."""
from typing import Optional
import pandas as pd

from data import(
    calculate_transmission_distances,
    get_years,
    set_transmission_tech_groups
    )

from utils import apply_dtypes

def set_transmission_losses(df_exist_corrected, df_planned_corrected,
                            centerpoints_dict, trn_param, SUBSEA_LINES):
    
    # Set base df including calculated transmission distances per pathway to set efficiencies.
    line_efficiency_dict = {}
    acdc_conv_efficiency_dict = {}

    # Get transmission parameters from config file.
    for tech, tech_params in trn_param.items():
        line_efficiency_dict[tech] = tech_params[2]
        acdc_conv_efficiency_dict[tech] = tech_params[3]
    
    # Assign technology groups to individual technologies.
    tech_group_dict = set_transmission_tech_groups(df_exist_corrected, 
                                                   df_planned_corrected, 
                                                   centerpoints_dict, 
                                                   trn_param, 
                                                   SUBSEA_LINES)
    
    # Set base df including calculated transmission distances per pathway.
    eff_df = calculate_transmission_distances(df_exist_corrected, df_planned_corrected, 
                                          centerpoints_dict)[['TECHNOLOGY', 'distance']]
    
    # Calculate technology and line specific transmission losses.
    eff_df['tech_group'] = eff_df['TECHNOLOGY'].map(tech_group_dict)
    
    for group in trn_param.keys():
        eff_df.loc[eff_df['tech_group'] == group, 
               'VALUE'] = round((eff_df['distance'] / 1000) * line_efficiency_dict[
                   group] + acdc_conv_efficiency_dict[group],3)
        
    # Convert losses into OAR format.
    eff_df['VALUE'] = 1 - (eff_df['VALUE'] / 100)
    
    return eff_df

def activity_transmission(df_iar_base, df_oar_base, df_eff,
                          start_year, end_year, region_name):

    # #### Downstream Activity Ratios
    
    # Build transmission system outputs
    
    df_iar_trn = df_oar_base.copy()
    
    # Change the technology name to PWRTRNXXXXX
    df_iar_trn["TECHNOLOGY"] = "PWRTRN" + df_iar_trn["FUEL"].str[3:8]
    # Make all modes of operation 1
    df_iar_trn["MODE_OF_OPERATION"] = 1
    # And remove all the duplicate entries
    df_iar_trn.drop_duplicates(keep="first", inplace=True)
    
    # Only keep technologies linked to the correct fuel
    df_iar_trn = df_iar_trn.loc[df_iar_trn['FUEL'].str.contains('ELC', na = False)]
    
    # OAR for transmission technologies is IAR, but the fuel is 02 instead of 01:
    df_oar_trn = df_iar_trn.copy()
    df_oar_trn["FUEL"] = df_oar_trn["FUEL"].str[0:8] + "02"
    
    # Build international transmission system from original input data, but for Line rather than Generator:    
    # Create MofO column:
    df_int_trn = df_eff[['TECHNOLOGY']].copy()
    df_int_trn["MODE_OF_OPERATION"] = 1
    df_int_trn = pd.concat([df_int_trn, df_int_trn.replace(int(1), int(2))])
    
    # Set the value (of either IAR or OAR) to 1 and add default columns.
    df_int_trn["VALUE"] = 1
    df_int_trn["REGION"] = region_name
    df_int_trn["YEAR"] = [get_years(start_year, end_year)] * len(df_int_trn)
    df_int_trn = df_int_trn.explode('YEAR').reset_index(drop = True)
    
    # Now create the input and output activity ratios
    df_int_trn_oar = df_int_trn.copy()
    df_int_trn_iar = df_int_trn.copy()
    
    # IAR Mode 1 is input from first country:
    df_int_trn_iar.loc[df_int_trn_iar["MODE_OF_OPERATION"] == 1, "FUEL"] = (
        "ELC" + df_int_trn_iar["TECHNOLOGY"].str[3:8] + "02"
    )
    # IAR Mode 2 is input from second country:
    df_int_trn_iar.loc[df_int_trn_iar["MODE_OF_OPERATION"] == 2, "FUEL"] = (
        "ELC" + df_int_trn_iar["TECHNOLOGY"].str[8:13] + "02"
    )

    # OAR Mode 2 is output to first country:
    df_int_trn_oar.loc[df_int_trn_oar["MODE_OF_OPERATION"] == 2, "FUEL"] = (
        "ELC" + df_int_trn_oar["TECHNOLOGY"].str[3:8] + "01"
    )
    # OAR Mode 1 is out to the second country:
    df_int_trn_oar.loc[df_int_trn_oar["MODE_OF_OPERATION"] == 1, "FUEL"] = (
        "ELC" + df_int_trn_oar["TECHNOLOGY"].str[8:13] + "01"
    )
    
    # Add line specific OAR values into matrix.
    df_int_trn_oar = df_int_trn_oar.drop(["VALUE"], axis=1)
    df_int_trn_oar = pd.merge(
        df_int_trn_oar, df_eff[['TECHNOLOGY', 'VALUE']], how="outer", on="TECHNOLOGY"
    )
    
    df_oar_trn_final = pd.concat(
            [
                df_oar_base,
                df_oar_trn,
                df_int_trn_oar,
            ]
        )
    
    df_iar_trn_final = pd.concat(
            [
                df_iar_base,
                df_iar_trn,
                df_int_trn_iar,
            ]
        )
    
    return df_iar_trn_final, df_oar_trn_final

def activity_transmission_limit(cross_border_trade, df_oar_trn_final):

    # Set cross-border trade to 0 if False
    if not cross_border_trade:
        df_crossborder_final = df_oar_trn_final[['REGION',
                                            'TECHNOLOGY'
                                            ]]
        df_crossborder_final = df_crossborder_final.drop_duplicates()
        df_crossborder_final = (df_crossborder_final
                                .loc[df_crossborder_final['TECHNOLOGY']
                                    .str.startswith('TRN')
                                    ]
                                )
        df_crossborder_final['VALUE'] = 0
    else:
        df_crossborder_final = pd.DataFrame(columns=['REGION', 
                                                    'TECHNOLOGY',
                                                    'VALUE'])

    df_crossborder_final = apply_dtypes(df_crossborder_final, 
                                        "TotalTechnologyModelPeriodActivityUpperLimit")
    
    return df_crossborder_final

def create_trn_dist_capacity_activity(*dfs: pd.DataFrame, value: Optional[float] = 31.536, 
                                      region: Optional[str] = "GLOBAL") -> pd.DataFrame:
    """Creates tranmission and distribution capacity to activity unit data
    
    Inputs are any number of dataframes with a TECHNOLOGY colums.
    """

    techs = []
    for df in dfs:
        temp = df[(df.TECHNOLOGY.str.startswith("TRN")) | (df.TECHNOLOGY.str.startswith("PWRTRN"))]
        techs += temp.TECHNOLOGY.to_list()
    
    data = []
    for tech in set(techs):
        data.append([region, tech, value])
        
    return pd.DataFrame(data, columns=["REGION", "TECHNOLOGY", "VALUE"])