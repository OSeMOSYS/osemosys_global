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