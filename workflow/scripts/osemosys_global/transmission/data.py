"""Functions to extract and format relevent data for tranmission."""
import pandas as pd

from sets import get_unique_technologies

def get_years(start: int, end: int) -> range:
    return range(start, end + 1)

def format_transmission_name(df):
    '''Formats PLEXOS transmission names into OSeMOSYS Global names.

    Args:
        :param df: Pandas DataFrame with columns 'From' and 'To' describing the 
               transmission from and to contries. ie. 
    
    Returns: 
        :param df: Same as df_in, except the 'From' and 'To' columns are replaced 
            with a single 'TECHNOLOGY' column holding OSeMOSYS Global 
            naming conventions 

    Example:
        df = pd.DataFrame(
            [[AF-COD, AF-COG, 0.001],
            [EU-AUT, EU-SVK, 0.004],
            [AS-LBN, AS-SYR, 0.006]], 
            columns = ['From', 'To', 'Losses']
        )
        pd.DataFrame(
            [[0.001,TRNCODXXCOGXX],
            [0.004,TRNAUTXXSVKXX],
            [0.006,TRNLBNXXSYRXX]] 
            columns = ['Losses','TECHNOLOGY'])'''

    # If from column has length 6 then it's the last three chars plus XX
    df.loc[df["From"].str.len() == 6, "From"] = (df["From"].str[3:6] + "XX")

    # If from column has length 9 then it's the 3:6 and 7:9 three chars plus XX
    df.loc[df["From"].str.len() == 9, "From"] = (
        df["From"].str[3:6] + df["From"].str[7:9])

    # If to column has length 6 then it's the last three chars plus XX
    df.loc[df["To"].str.len() == 6, "To"] = (df["To"].str[3:6] + "XX")

    # If to column has length 9 then it's the 3:6 and 7:9 three chars plus XX
    df.loc[df["To"].str.len() == 9, "To"] = (
        df["To"].str[3:6] + df["To"].str[7:9])

    # Combine From and To columns.
    df["TECHNOLOGY"] = ("TRN" + df["From"] + df["To"])

    # Drop to and from columns
    df = df.drop(["From", "To"], axis=1)

    return df

def format_gtd_existing(df):
    cols = {'pathway' : 'TECHNOLOGY',
            'from_region' : 'From',
            'to_region' : 'To', 
            'VALUE' : 'VALUE'
            }
    
    # Sets the capacity of the line to the average of max and max counter flows.
    df['VALUE'] = (df['max_flow'] + df['max_counter_flow']) / 2
    
    df = df[list(cols.keys())]
    
    df = df.rename(columns = cols)
    
    return df

def format_gtd_planned(df):
    cols = {'pathway' : 'TECHNOLOGY',
            'from_region' : 'From',
            'to_region' : 'To', 
            'year_planned' : 'YEAR',
            'VALUE' : 'VALUE'
            }
    
    # Sets the capacity of the line to the average of max and max counter flows.
    df['VALUE'] = (df['max_flow'] + df['max_counter_flow']) / 2
    
    df = df[list(cols.keys())]
    
    df = df.rename(columns = cols)
    
    return df

def set_gtd_mapping(df_exist, df_plan, region_mapping_dict, 
                         USA_TRN_BA_DICT_FROM, USA_TRN_BA_DICT_TO):
    
    # Create df with unique transmission technologies from GTD.
    tech_exist = get_unique_technologies(df_exist)
    tech_plan = get_unique_technologies(df_plan)
    tech_trn_list = list(set(tech_exist + tech_plan)) 
    df = pd.DataFrame(tech_trn_list, columns = ['TECHNOLOGY_gtd'])
    
    # Map OG regions to GTD technologies.
    df['From'], df['To'] = df["TECHNOLOGY_gtd"].str[3:8], df["TECHNOLOGY_gtd"].str[8:13]
    df['From'], df['To'] = df['From'].map(region_mapping_dict
                                          ), df['To'].map(region_mapping_dict)
    
    # For USA GTD regions that are larger than certain OG regions (MISO, PJM, SWPP),
    # assign custom OG regions as a best estimate for transmission landing points.
    df.loc[df['From'] == 'DUPLICATE', 
           'From'] = df['TECHNOLOGY_gtd'].map(USA_TRN_BA_DICT_FROM)
    
    df.loc[df['To'] == 'DUPLICATE', 
           'To'] = df['TECHNOLOGY_gtd'].map(USA_TRN_BA_DICT_TO)

    return df
