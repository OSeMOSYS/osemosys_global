"""Function to calaculate activity related to transmission."""

import pandas as pd

from constants import (
    region_name,
    start_year,
    end_year
)

from data_transmission import format_transmission_name

from utils import apply_dtypes

def activity_transmission(df_oar_base, df_pw_prop, df_trn_efficiencies):

    # #### Downstream Activity Ratios
    
    # Build transmission system outputs
    
    df_iar_trn = df_oar_base.copy()
    
    # Change the technology name to PWRTRNXXXXX
    df_iar_trn["TECHNOLOGY"] = "PWRTRN" + df_iar_trn["FUEL"].str[3:8]
    # Make all modes of operation 1
    df_iar_trn["MODE_OF_OPERATION"] = 1
    # And remove all the duplicate entries
    df_iar_trn.drop_duplicates(keep="first", inplace=True)
    
    # OAR for transmission technologies is IAR, but the fuel is 02 instead of 01:
    df_oar_trn = df_iar_trn.copy()
    df_oar_trn["FUEL"] = df_oar_trn["FUEL"].str[0:8] + "02"
    
    # Build international transmission system from original input data, but for Line rather than Generator:
    int_trn_cols = ["child_class", "child_object", "property", "value"]
    df_int_trn = df_pw_prop[int_trn_cols]
    df_int_trn = df_int_trn[df_int_trn["child_class"] == "Line"]
    
    # For IAR and OAR we can drop the value:
    df_int_trn = df_int_trn.drop(["child_class", "value"], axis=1)
    
    # Create MofO column based on property:
    df_int_trn["MODE_OF_OPERATION"] = 1
    df_int_trn.loc[df_int_trn["property"] == "Min Flow", "MODE_OF_OPERATION"] = 2
    
    # Use the child_object column to build the technology names:
    df_int_trn["codes"] = df_int_trn["child_object"].str.split(pat="-")
    
    # If there are only two locations, then the node is XX
    df_int_trn.loc[df_int_trn["codes"].str.len() == 2, "TECHNOLOGY"] = (
        "TRN" + df_int_trn["codes"].str[0] + "XX" + df_int_trn["codes"].str[1] + "XX"
    )
    # If there are four locations, the node is already included
    df_int_trn.loc[df_int_trn["codes"].str.len() == 4, "TECHNOLOGY"] = (
        "TRN"
        + df_int_trn["codes"].str[0]
        + df_int_trn["codes"].str[1]
        + df_int_trn["codes"].str[2]
        + df_int_trn["codes"].str[3]
    )
    # If there are three items, and the last item is two characters, then the second item is an XX:
    df_int_trn.loc[
        (df_int_trn["codes"].str.len() == 3) & (df_int_trn["codes"].str[2].str.len() == 2),
        "TECHNOLOGY",
    ] = (
        "TRN"
        + df_int_trn["codes"].str[0]
        + "XX"
        + df_int_trn["codes"].str[1]
        + df_int_trn["codes"].str[2]
    )
    # If there are three items, and the last item is three characters, then the last item is an XX:
    df_int_trn.loc[
        (df_int_trn["codes"].str.len() == 3) & (df_int_trn["codes"].str[2].str.len() == 3),
        "TECHNOLOGY",
    ] = (
        "TRN"
        + df_int_trn["codes"].str[0]
        + df_int_trn["codes"].str[1]
        + df_int_trn["codes"].str[2]
        + "XX"
    )
    
    # Set the value (of either IAR or OAR) to 1
    df_int_trn["VALUE"] = 1
    df_int_trn["REGION"] = region_name
    
    df_int_trn = df_int_trn.drop(["property", "child_object", "codes"], axis=1)
    df_int_trn["YEAR"] = start_year
    
    # add in future years
    df_int_trn_new = df_int_trn.copy()
    df_int_trn_new["YEAR"] = [range(start_year + 1, end_year + 1)] * len(df_int_trn_new)
    df_int_trn_new = df_int_trn_new.explode(column="YEAR")
    
    df_int_trn = pd.concat([df_int_trn, df_int_trn_new]).reset_index(drop=True)
    
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
    
    # Drop unneeded columns
    df_trn_efficiencies = df_trn_efficiencies.drop(
        [
            "Interface",
            "KM distance",
            "HVAC/HVDC/Subsea",
            "Build Cost ($2010 in $000)",
            "Build cost + FOM (3.5% of Capex per year)",
            "Unnamed: 8",
            "Interface is per MW",
            "Unnamed: 10",
            "Unnamed: 11",
            "Unnamed: 12",
            "Subsea lines",
            "Unnamed: 14"
        ],
        axis=1,
    )
    
    # Drop NaN values
    df_trn_efficiencies = df_trn_efficiencies.dropna(subset=["From"])
    
    # Create tech column from To and From Codes:
    df_trn_efficiencies = format_transmission_name(df_trn_efficiencies)
    
    # Rename column 'VALUES'
    df_trn_efficiencies = df_trn_efficiencies.rename(columns={"Losses": "VALUE"})
    
    # And adjust OAR values to be output amounts vs. losses:
    df_trn_efficiencies['VALUE'] = 1.0 - df_trn_efficiencies['VALUE']
    
    # and add values into OAR matrix
    df_int_trn_oar = df_int_trn_oar.drop(["VALUE"], axis=1)
    df_int_trn_oar = pd.merge(
        df_int_trn_oar, df_trn_efficiencies, how="outer", on="TECHNOLOGY"
    )
    
    return df_iar_trn, df_oar_trn, df_int_trn_oar, df_int_trn_iar

def activity_transmission_limit(cross_border_trade, df_oar_final):

    # Set cross-border trade to 0 if False
    if not cross_border_trade:
        df_crossborder_final = df_oar_final[['REGION',
                                            'TECHNOLOGY'
                                            ]]
        df_crossborder_final = df_crossborder_final.drop_duplicates()
        df_crossborder_final = (df_crossborder_final
                                .loc[df_crossborder_final['TECHNOLOGY']
                                    .str.startswith('TRN')
                                    ]
                                )
        df_crossborder_final['VALUE'] = 0
    else:
        df_crossborder_final = pd.DataFrame(columns=['REGION', 
                                                    'TECHNOLOGY',
                                                    'VALUE'])

    df_crossborder_final = apply_dtypes(df_crossborder_final, 
                                        "TotalTechnologyModelPeriodActivityUpperLimit")
    
    return df_crossborder_final