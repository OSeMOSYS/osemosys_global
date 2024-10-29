"""Module for reading in data sources"""

import pandas as pd

def import_op_life(f: str) -> pd.DataFrame:
    """Imports default operational life data.
    
    operational_life.csv
    """
    return pd.read_csv(f)

def import_iar_base(f: str) -> pd.DataFrame:
    """Imports InputActivityRatio.csv as output from the Transmission rule.
    
    InputActivityRatio.csv
    """
    return pd.read_csv(f)

def import_oar_base(f: str) -> pd.DataFrame:
    """Imports OutputActivityRatio.csv as output from the Transmission rule.
    
    OutputActivityRatio.csv
    """
    return pd.read_csv(f)

def import_capact_base(f: str) -> pd.DataFrame:
    """Imports CapacityToActivityUnit.csv as output from the Transmission rule.
    
    CapacityToActivityUnit.csv
    """
    return pd.read_csv(f)

def import_fix_cost_base(f: str) -> pd.DataFrame:
    """Imports FixedCost.csv as output from the Transmission rule.
    
    FixedCost.csv
    """
    return pd.read_csv(f)

def import_var_cost_base(f: str) -> pd.DataFrame:
    """Imports VariableCost.csv as output from the Transmission rule.
    
    FixedCost.csv
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