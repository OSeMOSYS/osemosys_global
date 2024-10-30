"""Module for reading in data sources"""

import pandas as pd

def import_gtd_existing(f: str) -> pd.DataFrame:
    """Imports existing transmission capacity data from the Global
    Transmission Database (GTD)."""
    return pd.read_csv(f, encoding = 'latin1')

def import_gtd_planned(f: str) -> pd.DataFrame:
    """Imports planned transmission capacity data from the Global
    Transmission Database (GTD)."""
    return pd.read_csv(f, encoding = 'latin1')

def import_gtd_mapping(f: str) -> pd.DataFrame:
    """Imports the spatial mapping between OG and the Global
    Transmission Database (GTD)."""
    return pd.read_csv(f)

def import_centerpoints(f: str) -> pd.DataFrame:
    """Imports the centerpoints as used to calculate transmission
    distances for transmission lines."""
    return pd.read_csv(f)

def import_transmission_build_rates(f: str) -> pd.DataFrame:
    """Imports transmission pathway specific user defined max build rates."""
    return pd.read_csv(f)

def import_op_life(f: str) -> pd.DataFrame:
    """Imports default operational life data.
    
    operational_life.csv
    """
    return pd.read_csv(f)

def import_iar_base(f: str) -> pd.DataFrame:
    """Imports InputActivityRatio.csv as output from the Powerplant rule.
    
    InputActivityRatio.csv
    """
    return pd.read_csv(f)

def import_oar_base(f: str) -> pd.DataFrame:
    """Imports OutputActivityRatio.csv as output from the Powerplant rule.
    
    OutputActivityRatio.csv
    """
    return pd.read_csv(f)

def import_capact_base(f: str) -> pd.DataFrame:
    """Imports CapacityToActivityUnit.csv as output from the Powerplant rule.
    
    CapacityToActivityUnit.csv
    """
    return pd.read_csv(f)

def import_cap_cost_base(f: str) -> pd.DataFrame:
    """Imports CapitalCost.csv as output from the Powerplant rule.
    
    CapitalCost.csv
    """
    return pd.read_csv(f)

def import_fix_cost_base(f: str) -> pd.DataFrame:
    """Imports FixedCost.csv as output from the Powerplant rule.
    
    FixedCost.csv
    """
    return pd.read_csv(f)

def import_op_life_base(f: str) -> pd.DataFrame:
    """Imports OperationalLife.csv as output from the Powerplant rule.
    
    OperationalLife.csv
    """
    return pd.read_csv(f)

def import_max_cap_invest_base(f: str) -> pd.DataFrame:
    """Imports TotalAnnualMaxCapacityInvestment.csv as output from the Powerplant rule.
    
    TotalAnnualMaxCapacityInvestment.csv
    """
    return pd.read_csv(f)

def import_min_cap_invest_base(f: str) -> pd.DataFrame:
    """Imports TotalAnnualMinCapacityInvestment.csv as output from the Powerplant rule.
    
    TotalAnnualMinCapacityInvestment.csv
    """
    return pd.read_csv(f)

def import_res_cap_base(f: str) -> pd.DataFrame:
    """Imports ResidualCapacity.csv as output from the Powerplant rule.
    
    ResidualCapacity.csv
    """
    return pd.read_csv(f)

def import_set_base(f: str) -> pd.DataFrame:
    """Imports a set csv"""
    return pd.read_csv(f)