"""Functions to calaculate activity related to transmission."""
from typing import Optional
import pandas as pd

from data import(
    calculate_transmission_distances,
    get_years,
    set_transmission_tech_groups
    )

from utils import apply_dtypes

# %% 

def activity_storage(storage_set, df_iar_base, df_oar_base, df_eff,
                     start_year, end_year, region_name):

    # Add InputActivityRatio and OutputActivityRatio
    # InputActivityRatio
    df_storage_iar = pd.DataFrame(
        list(itertools.product([region_name], storage_techs, years, [1])),
        columns=["REGION", "TECHNOLOGY", "YEAR", "MODE_OF_OPERATION"],
    )
    df_storage_iar["VALUE"] = 1
    df_storage_iar["FUEL"] = "ELC" + df_storage_iar["TECHNOLOGY"].str[6:11] + "01"
    df_storage_iar = df_storage_iar[
        ["REGION", "TECHNOLOGY", "FUEL", "MODE_OF_OPERATION", "YEAR", "VALUE"]
    ]
    
    wait_time = 0
    while not os.path.exists(os.path.join(output_data_dir, "InputActivityRatio.csv")):
        time.sleep(5)
        wait_time += 1
        if wait_time > 20:
            break
    df_iar = pd.read_csv(os.path.join(output_data_dir, "InputActivityRatio.csv"))
    df_iar = pd.concat([df_iar, df_storage_iar])
    df_iar.to_csv(os.path.join(output_data_dir, "InputActivityRatio.csv"), index=None)
    time.sleep(20)
    
    # OutputActivityRatio
    df_storage_oar = pd.DataFrame(
        list(itertools.product([region_name], storage_techs, years, [2])),
        columns=["REGION", "TECHNOLOGY", "YEAR", "MODE_OF_OPERATION"],
    )
    df_storage_oar["VALUE"] = 1
    df_storage_oar["FUEL"] = "ELC" + df_storage_oar["TECHNOLOGY"].str[6:11] + "01"
    df_storage_oar = df_storage_oar[
        ["REGION", "TECHNOLOGY", "FUEL", "MODE_OF_OPERATION", "YEAR", "VALUE"]
    ]
    
    wait_time = 0
    while not os.path.exists(os.path.join(output_data_dir, "OutputActivityRatio.csv")):
        time.sleep(5)
        wait_time += 1
        if wait_time > 20:
            break
    df_oar = pd.read_csv(os.path.join(output_data_dir, "OutputActivityRatio.csv"))
    df_oar = pd.concat([df_oar, df_storage_oar])
    df_oar.to_csv(os.path.join(output_data_dir, "OutputActivityRatio.csv"), index=None)
    time.sleep(20)


# %%

def activity_transmission(df_iar_base, df_oar_base, df_eff,
                          start_year, end_year, region_name):

    # #### Downstream Activity Ratios
    
    # Build transmission system outputs
    
    df_iar_trn = df_oar_base.copy()
    
    # Change the technology name to PWRTRNXXXXX
    df_iar_trn["TECHNOLOGY"] = "PWRTRN" + df_iar_trn["FUEL"].str[3:8]
    # Make all modes of operation 1
    df_iar_trn["MODE_OF_OPERATION"] = 1
    # And remove all the duplicate entries
    df_iar_trn.drop_duplicates(keep="first", inplace=True)
    
    # Only keep technologies linked to the correct fuel
    df_iar_trn = df_iar_trn.loc[df_iar_trn['FUEL'].str.contains('ELC', na = False)]
    
    # OAR for transmission technologies is IAR, but the fuel is 02 instead of 01:
    df_oar_trn = df_iar_trn.copy()
    df_oar_trn["FUEL"] = df_oar_trn["FUEL"].str[0:8] + "02"
    
    # Build international transmission system from original input data, but for Line rather than Generator:    
    # Create MofO column:
    df_int_trn = df_eff[['TECHNOLOGY']].copy()
    df_int_trn["MODE_OF_OPERATION"] = 1
    df_int_trn = pd.concat([df_int_trn, df_int_trn.replace(int(1), int(2))])
    
    # Set the value (of either IAR or OAR) to 1 and add default columns.
    df_int_trn["VALUE"] = 1
    df_int_trn["REGION"] = region_name
    df_int_trn["YEAR"] = [get_years(start_year, end_year)] * len(df_int_trn)
    df_int_trn = df_int_trn.explode('YEAR').reset_index(drop = True)
    
    # Now create the input and output activity ratios
    df_int_trn_oar = df_int_trn.copy()
    df_int_trn_iar = df_int_trn.copy()
    
    # IAR Mode 1 is input from first country:
    df_int_trn_iar.loc[df_int_trn_iar["MODE_OF_OPERATION"] == 1, "FUEL"] = (
        "ELC" + df_int_trn_iar["TECHNOLOGY"].str[3:8] + "02"
    )
    # IAR Mode 2 is input from second country:
    df_int_trn_iar.loc[df_int_trn_iar["MODE_OF_OPERATION"] == 2, "FUEL"] = (
        "ELC" + df_int_trn_iar["TECHNOLOGY"].str[8:13] + "02"
    )

    # OAR Mode 2 is output to first country:
    df_int_trn_oar.loc[df_int_trn_oar["MODE_OF_OPERATION"] == 2, "FUEL"] = (
        "ELC" + df_int_trn_oar["TECHNOLOGY"].str[3:8] + "01"
    )
    # OAR Mode 1 is out to the second country:
    df_int_trn_oar.loc[df_int_trn_oar["MODE_OF_OPERATION"] == 1, "FUEL"] = (
        "ELC" + df_int_trn_oar["TECHNOLOGY"].str[8:13] + "01"
    )
    
    # Add line specific OAR values into matrix.
    df_int_trn_oar = df_int_trn_oar.drop(["VALUE"], axis=1)
    df_int_trn_oar = pd.merge(
        df_int_trn_oar, df_eff[['TECHNOLOGY', 'VALUE']], how="outer", on="TECHNOLOGY"
    )
    
    df_oar_trn_final = pd.concat(
            [
                df_oar_base,
                df_oar_trn,
                df_int_trn_oar,
            ]
        )
    
    df_iar_trn_final = pd.concat(
            [
                df_iar_base,
                df_iar_trn,
                df_int_trn_iar,
            ]
        )
    
    return df_iar_trn_final, df_oar_trn_final

def create_trn_dist_capacity_activity(*dfs: pd.DataFrame, value: Optional[float] = 31.536, 
                                      region: Optional[str] = "GLOBAL") -> pd.DataFrame:
    """Creates tranmission and distribution capacity to activity unit data
    
    Inputs are any number of dataframes with a TECHNOLOGY colums.
    """

    techs = []
    for df in dfs:
        temp = df[(df.TECHNOLOGY.str.startswith("TRN")) | (df.TECHNOLOGY.str.startswith("PWRTRN"))]
        techs += temp.TECHNOLOGY.to_list()
    
    data = []
    for tech in set(techs):
        data.append([region, tech, value])
        
    return pd.DataFrame(data, columns=["REGION", "TECHNOLOGY", "VALUE"])