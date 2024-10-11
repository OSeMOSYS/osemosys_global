"""Functions to extract and format relevent data for tranmission."""
import pandas as pd
import geopy.distance

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

def correct_gtd_line_mapping(df_exist, df_plan, region_mapping_dict, 
                             CUSTOM_TRN_BA_DICT_FROM, CUSTOM_TRN_BA_DICT_TO):
    
    # Create df with unique transmission technologies from GTD.
    tech_exist = get_unique_technologies(df_exist)
    tech_plan = get_unique_technologies(df_plan)
    tech_trn_list = list(set(tech_exist + tech_plan)) 
    df = pd.DataFrame(tech_trn_list, columns = ['TECHNOLOGY'])
    
    # Map OG regions to GTD technologies.
    df['From'], df['To'] = df['TECHNOLOGY'].str[3:8], df['TECHNOLOGY'].str[8:13]
    df['From'], df['To'] = df['From'].map(region_mapping_dict
                                          ), df['To'].map(region_mapping_dict)
    
    # For USA GTD regions that are larger than certain OG regions (MISO, PJM, SWPP),
    # as well as for the Nei Mongol region in China, assign custom OG regions as a 
    # best estimate for transmission landing points.
    df.loc[df['From'] == 'DUPLICATE', 
           'From'] = df['TECHNOLOGY'].map(CUSTOM_TRN_BA_DICT_FROM)
    
    df.loc[df['To'] == 'DUPLICATE', 
           'To'] = df['TECHNOLOGY'].map(CUSTOM_TRN_BA_DICT_TO)
    
    return df
    
def calculate_transmission_distances(df_exist, df_plan, region_mapping_dict, 
                                     CUSTOM_TRN_BA_DICT_FROM, CUSTOM_TRN_BA_DICT_TO, 
                                     centerpoints_dict):
    
    df = correct_gtd_line_mapping(df_exist, df_plan, region_mapping_dict,
                                  CUSTOM_TRN_BA_DICT_FROM, CUSTOM_TRN_BA_DICT_TO)
    
    for row in centerpoints_dict:
      df.loc[df['From'] == row['region'], 'From_lat'] = row['lat']
      df.loc[df['From'] == row['region'], 'From_long'] = row['long']
      
      df.loc[df['To'] == row['region'], 'To_lat'] = row['lat']
      df.loc[df['To'] == row['region'], 'To_long'] = row['long']
    
    df['distance'] = df.apply(lambda x: geopy.distance.geodesic((x['From_lat'], x['From_long']), 
                                                          (x['To_lat'], x['To_long'])), axis=1)
    
    df['distance'] = df['distance'].astype(str).str.replace(' km', '').astype(float).round(0)
    
    return df