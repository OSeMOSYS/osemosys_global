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

def import_res_limit(f: str) -> pd.DataFrame:
    """Imports the PLEXOS-World MESSAGix soft link model file for 
    renewable resource limits.
    
    PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx
    """
    return pd.read_excel(f, sheet_name = 'Properties')

def import_fuel_limit(f: str) -> pd.DataFrame:
    """Imports user defined fuel limits for powerplant technologies.
    
    fuel_limits.csv
    """
    return pd.read_csv(f)

def import_build_rates(f: str) -> pd.DataFrame:
    """Imports user defined build rates for powerplant technologies.
    
    powerplant_build_rates.csv
    """
    return pd.read_csv(f)

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

def import_custom_res_potentials(f: str) -> pd.DataFrame:
    """Imports renewable potentials for custom nodes.
    
    custom_nodes\RE_potentials.csv
    """
    return pd.read_csv(f)

def import_cmo_forecasts(f: str) -> pd.DataFrame:
    """Imports CMO forecasts.
    
    CMO-October-2024-Forecasts.csv
    """
    return pd.read_excel(f, header=1, skiprows=(
        [i for i in range(1, 3)] + [j for j in range(4, 25)]), nrows=5)

def import_fuel_prices(f: str) -> pd.DataFrame:
    """Imports international and country specific fuel prices.
    
    fuel_prices.csv
    """
    return pd.read_csv(f)

def import_specified_annual_demand(f: str) -> pd.DataFrame:
    """Imports SpecifiedAnnualDemand.csv as output from the demand_projections rule.
    
    SpecifiedAnnualDemand.csv
    """
    return pd.read_csv(f)