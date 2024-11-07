"""Functions to set the technology to/from storage files."""
import pandas as pd
import itertools

def set_base_df(ar, region_name):
    
    data = pd.DataFrame(
        columns=["REGION", "TECHNOLOGY", "STORAGE", "MODE_OF_OPERATION"])
    
    for sto in ar['VALUE']:
        df = pd.DataFrame(
            list(
                itertools.product(
                    [region_name],
                    ["PWR" + sto],
                    [sto],
                    [1, 2],
                )
            ),
            columns=["REGION", "TECHNOLOGY", "STORAGE", "MODE_OF_OPERATION"],
        )
        
        data = pd.concat([data, df])
        
    return data

def set_technology_to_storage(ar, region_name):
    data = set_base_df(ar, region_name)
    data.loc[data["MODE_OF_OPERATION"] == 1, "VALUE"] = 1.0
    data.loc[data["MODE_OF_OPERATION"] == 2, "VALUE"] = 0.0
    data["VALUE"] = data["VALUE"].astype(float)
    
    return data
    
def set_technology_from_storage(ar, region_name):
    data = set_base_df(ar, region_name)
    data.loc[data["MODE_OF_OPERATION"] == 1, "VALUE"] = 0.0
    data.loc[data["MODE_OF_OPERATION"] == 2, "VALUE"] = 1.0
    data["VALUE"] = data["VALUE"].astype(float)
    
    return data