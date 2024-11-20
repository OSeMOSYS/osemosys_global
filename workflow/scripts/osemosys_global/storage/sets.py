"""Function to create sets."""

import pandas as pd

def set_unique_storage_technologies(ar: pd.DataFrame, sto_techs: list) -> pd.DataFrame():
    ar = [
        x for x in ar["VALUE"].unique().astype(str) if 
        x.startswith("ELC") if x.endswith("01")
    ]
    
    data = []
    
    for tech in sto_techs:
        for val in ar:
            data.append(val.replace('ELC', tech))
    
    return pd.Series(data, name='VALUE').to_frame()

def set_unique_technologies(ar: pd.DataFrame) -> pd.DataFrame():
    data = ar.copy()
    data['VALUE'] = 'PWR' + data['VALUE']
    
    return data