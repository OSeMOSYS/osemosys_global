"""Functions to calaculate activity related to transmission."""
from typing import Optional
import pandas as pd
import itertools

from data import get_years

def activity_storage(storage_set, df_iar_base, df_oar_base, storage_param,
                     start_year, end_year, region_name):
    
    efficiency_dict = {}
    
    # Set baseline technology efficiencies as defined in the config file.
    for tech, tech_params in storage_param.items():
        efficiency_dict[tech] = tech_params[1]
    
    years = [get_years(start_year, end_year)]

    # InputActivityRatio
    df_storage_iar = pd.DataFrame(
        list(itertools.product([region_name], storage_set['VALUE'].unique(), years, [1])),
        columns=["REGION", "TECHNOLOGY", "YEAR", "MODE_OF_OPERATION"],
    )
    df_storage_iar = df_storage_iar.explode('YEAR').reset_index(drop = True)
    df_storage_iar["VALUE"] = 1
    df_storage_iar["TECHNOLOGY"] = "PWR" + df_storage_iar["TECHNOLOGY"]  
    df_storage_iar["FUEL"] = "ELC" + df_storage_iar["TECHNOLOGY"].str[3:8] + "01"
    df_storage_iar = df_storage_iar[
        ["REGION", "TECHNOLOGY", "FUEL", "MODE_OF_OPERATION", "YEAR", "VALUE"]
    ]
    
    df_iar = pd.concat([df_iar_base, df_storage_iar])

    # OutputActivityRatio
    df_storage_oar = pd.DataFrame(
        list(itertools.product([region_name], storage_set['VALUE'].unique(), years, [2])),
        columns=["REGION", "TECHNOLOGY", "YEAR", "MODE_OF_OPERATION"],
    )
    df_storage_oar = df_storage_oar.explode('YEAR').reset_index(drop = True)
    
    for each_tech in storage_param.keys():
        df_storage_oar.loc[df_storage_oar['TECHNOLOGY'].str.contains(each_tech),
                           'VALUE'] = round(1 / (efficiency_dict[each_tech] / 100), 3)

    df_storage_oar["TECHNOLOGY"] = "PWR" + df_storage_oar["TECHNOLOGY"] 
    df_storage_oar["FUEL"] = "ELC" + df_storage_oar["TECHNOLOGY"].str[3:8] + "01"
    df_storage_oar = df_storage_oar[
        ["REGION", "TECHNOLOGY", "FUEL", "MODE_OF_OPERATION", "YEAR", "VALUE"]
    ]
    
    df_oar = pd.concat([df_oar_base, df_storage_oar])
    
    return df_iar, df_oar

def create_storage_capacity_activity(storage_set: pd.DataFrame, df_capact_base: pd.DataFrame, 
                                     value: Optional[float] = 31.536, 
                                     region: Optional[str] = "GLOBAL") -> pd.DataFrame:
    """Creates and concats storage to activity unit data
    """
    
    data = []
    for tech in storage_set['VALUE'].unique():
        data.append([region, 'PWR' + tech, value])
        
    data = pd.DataFrame(data, columns=["REGION", "TECHNOLOGY", "VALUE"])
        
    data = pd.concat([df_capact_base, data])
        
    return data