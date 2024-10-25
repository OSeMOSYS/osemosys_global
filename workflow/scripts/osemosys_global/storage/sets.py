"""Function to create sets."""

import pandas as pd

def set_unique_storage_technologies(ar: pd.DataFrame, sto_techs: list) -> pd.DataFrame():
    ar = [
        x for x in ar["VALUE"].unique() if x.startswith("ELC") if x.endswith("01")
    ]
    
    data = []
    
    for tech in sto_techs:
        for val in ar:
            data.append(val.replace('ELC', tech))
    
    return pd.Series(data, name='VALUE').to_frame()

def set_unique_storage_technologies_custom_nodes(nodes: list, sto_techs: list) -> pd.DataFrame():
    data = []
    
    for tech in sto_techs:
        for val in nodes:
            data.append(tech + val + '01')
    
    return pd.Series(data, name="VALUE").to_frame()

def set_unique_technologies(ar: pd.DataFrame) -> pd.DataFrame():
    data = ar.copy()
    data['VALUE'] = 'PWR' + data['VALUE']
    
    return data

#def create_set_from_iterators(*args: list[str] | set[str]) -> pd.DataFrame:
#    """Creates a set formatted dataframe from an arbitrary number of arguments"""
    
#    data = []
    
  #  for arg in args:
  #      data += arg
 #   
  #  data = list(set(data))
    
  #  return pd.Series(data, name="VALUE").to_frame().sort_values(['VALUE'])