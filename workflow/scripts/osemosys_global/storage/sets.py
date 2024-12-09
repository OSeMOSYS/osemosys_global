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

def set_storage_fuels(ar: pd.DataFrame) -> pd.DataFrame():
    data = ar.copy()
    data = data.loc[(data['VALUE'].str.startswith('ELC')) & 
                    (data['VALUE'].str.endswith('01'))]
    
    data['VALUE'] = data['VALUE'].str.replace('01', '03')
    data = pd.concat([ar, data])
    
    return data