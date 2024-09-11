"""Functions to set Spatial Mapping from PLEXOS"""

import pandas as pd

def _get_spatial_mapping(plexos: pd.DataFrame) -> pd.DataFrame:
    """Gets general spatial mapping data structure"""
    
    df = plexos.copy()
    df = df[df.collection.str.contains("Region")].copy()
    df["Country"] = df.parent_object.str.split("-", expand=True)[1]
    return df

def get_spatial_mapping_country(plexos: pd.DataFrame) -> pd.DataFrame:
    """Gets spatial mapping by country"""
    
    df = _get_spatial_mapping(plexos).set_index("Country")
    df = df.loc[~df.index.duplicated(keep="first")]
    return df
    
def get_spatial_mapping_node(plexos: pd.DataFrame) -> pd.DataFrame:
    """Gets spatial mapping by node"""
    
    df = _get_spatial_mapping(plexos).set_index("parent_object")
    df = df.loc[~df.index.duplicated(keep="first")]
    return df
    