"""Module for reading in data sources"""

import pandas as pd

def import_plexos_2015(f: str, metric: str) -> dict[str, pd.DataFrame]:
    """Imports PLEXOS-World 2015 model file.
    
    PLEXOS_World_2015_Gold_V1.1.xlsx
    """
    if metric.lower() == "memb":
        sheet_name = "Memberships"
    elif metric.lower() == "prop":
        sheet_name = "Properties"
    else:
        raise NotImplementedError

    return pd.read_excel(f, sheet_name=sheet_name)

def import_weo_regions(f: str) -> pd.DataFrame:
    """Imports WEO to OG region mapping.
    
    weo_region_mapping.csv
    """
    return pd.read_csv(f)

def import_weo_costs(f: str) -> pd.DataFrame:
    """Imports WEO cost data for powerplant technologies.
    
    weo_2020_powerplant_costs.csv
    """
    return pd.read_csv(f)

def import_op_life(f: str) -> pd.DataFrame:
    """Imports default operational life data.
    
    operational_life.csv
    """
    return pd.read_csv(f)

def import_naming_convention_tech(f: str) -> pd.DataFrame:
    """Imports naming convention for powerplant technologies.
    
    naming_convention_tech.csv
    """
    return pd.read_csv(f)

def import_line_data(f: str, metric: str) -> dict[str, pd.DataFrame]:
    """Imports transmission data from PLEXOS-World.
    
    Costs Line expansion.xlsx
    """
    if metric.lower() == "interface":
        sheet_name = "Interface"
    elif metric.lower() == "lines":
        sheet_name = "Lines"
    else:
        raise NotImplementedError

    return pd.read_excel(f, sheet_name=sheet_name)

def import_afs(f: str) -> pd.DataFrame:
    """Imports default availability factors.
    
    availability_factors.csv
    """
    return pd.read_csv(f)

def import_custom_res_cap(f: str) -> pd.DataFrame:
    """Imports residual capacities for custom nodes.
    
    custom_nodes\residual_capacity.csv
    """
    return pd.read_csv(f)

def import_tech_set(f: str) -> pd.DataFrame:
    """Imports list of technologies.
    
    TECHNOLOGY.csv
    """
    return pd.read_csv(f)

def import_fuel_set(f: str) -> pd.DataFrame:
    """Imports list of fuels.
    
    FUEL.csv
    """
    return pd.read_csv(f)