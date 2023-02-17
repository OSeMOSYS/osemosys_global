"""Module for utility plotting functions."""

import pandas as pd
import os 
from osemosys_global.OPG_configuration import ConfigFile, ConfigPaths
from typing import Dict
import itertools
from osemosys_global.visualisation.constants import DAYS_PER_MONTH, MONTH_NAMES

def get_color_codes() -> Dict:
    """Get color naming dictionary.
    
    Return:
        Dictionary with tech and color name map
    """
    config_paths = ConfigPaths()
    input_data_dir = config_paths.input_data_dir
    name_colour_codes = pd.read_csv(os.path.join(input_data_dir,
                                                'color_codes.csv'
                                                ),
                                   encoding='latin-1')

    # Get colour mapping dictionary
    color_dict = dict([(n, c) for n, c
                   in zip(name_colour_codes.tech_id,
                          name_colour_codes.colour)])
    return color_dict

def get_regions(
    data: dict[str, pd.DataFrame],
    countries_only: bool = False,
) -> list:
    """Gets all region codes
    
    Arguments:
        data: dict[str, pd.DataFrame]
            Input datastore 
        countries_only: bool
            Get only country if true, else get 5 letter region code
    
    Returns
        List of country/region codes
    """
    region_data = region_filter(data["TECHNOLOGY"])
    if countries_only:
        return list(set([x[0:3] for x in region_data]))
    else:
        return region_data
    
def region_filter(df:pd.DataFrame) -> list:
    """Filters out unique region codes from TECHNOLOGY codes.
    
    Arguments
        df:pd.DataFrame
            Dataframe with a TECHNOLOGY columns
    
    Returns
        list
            All regions emcompased in the TECHNOLOGY column 
    """
    techs = df.loc[
        (df["VALUE"].str.startswith("PWR")) &
        ~(df["VALUE"].str.startswith("PWRTRN"))]
    regions = techs["VALUE"].str[-7:-2]
    return regions.drop_duplicates().to_list()

def powerplant_filter(df: pd.DataFrame, country:str = None) -> pd.DataFrame:
    """Filters dataframe by powerplant generator (PWR)
    
    Arguments: 
        df: pd.DataFrame
            Input result dataframe  
        country: str 
            If a country provided, get data at a country level, else get data
            at a system level
    """

    filtered_df = df[~df.TECHNOLOGY.str.contains('TRN')]
    filtered_df = filtered_df.loc[filtered_df.TECHNOLOGY.str[0:3] == 'PWR']
    filtered_df['TYPE'] = filtered_df.TECHNOLOGY.str[3:6]
    filtered_df['COUNTRY'] = filtered_df.TECHNOLOGY.str[6:9]

    if country:
        filtered_df = filtered_df.loc[filtered_df['COUNTRY'] == country]
        filtered_df['LABEL'] = filtered_df['COUNTRY'] + '-' + filtered_df['TYPE']
    else:
        filtered_df['LABEL'] = filtered_df['TYPE']
    
    filtered_df.drop(['TECHNOLOGY', 'TYPE', 'COUNTRY'],
            axis=1,
            inplace=True)
    return filtered_df

def transform_ts(data:dict[str, pd.DataFrame], df:pd.DataFrame) -> pd.DataFrame:
    """Adds month, hour, year columns to timesliced data. 
    
    Arguments:
        data: dict[str, pd.DataFrame]
            Input datastore 
        df: pd.DataFrame
            Dataframe to add ts data to 
        
    Returns: 
        pd.DataFrame 
    """

    generation = list(data["TECHNOLOGY"]["VALUE"].unique())

    config = ConfigFile('config')
    seasons_raw = config.get('seasons')
    seasonsData = []

    for s, months in seasons_raw.items():
        for month in months:
            seasonsData.append([month, s]) 
    seasons_df = pd.DataFrame(seasonsData, 
                              columns=['month', 'season'])
    seasons_df = seasons_df.sort_values(by=['month']).reset_index(drop=True)
    dayparts_raw = config.get('dayparts')
    daypartData = []
    for dp, hr in dayparts_raw.items():
        daypartData.append([dp, hr[0], hr[1]])
    dayparts_df = pd.DataFrame(daypartData,
                               columns=['daypart', 'start_hour', 'end_hour'])
    timeshift = config.get('timeshift')
    dayparts_df['start_hour'] = dayparts_df['start_hour'].map(lambda x: apply_timeshift(x, timeshift))
    dayparts_df['end_hour'] = dayparts_df['end_hour'].map(lambda x: apply_timeshift(x, timeshift))

    month_names = MONTH_NAMES
    days_per_month = DAYS_PER_MONTH

    seasons_df['month_name'] = seasons_df['month'].map(month_names)
    seasons_df['days'] = seasons_df['month_name'].map(days_per_month)
    seasons_df_grouped = seasons_df.groupby(['season'],
                                            as_index=False)['days'].sum()
    days_dict = dict(zip(list(seasons_df_grouped['season']),
                         list(seasons_df_grouped['days'])
                         )
                     )
    seasons_df['days'] = seasons_df['season'].map(days_dict)

    years = config.get_years()

    seasons_dict = dict(zip(list(seasons_df['month']),
                            list(seasons_df['season'])
                            )
                        )

    dayparts_dict = {i: [j, k]
                     for i, j, k
                     in zip(list(dayparts_df['daypart']),
                            list(dayparts_df['start_hour']),
                            list(dayparts_df['end_hour'])
                            )
                     }

    hours_dict = {i: abs(k-j)
                  for i, j, k
                  in zip(list(dayparts_df['daypart']),
                         list(dayparts_df['start_hour']),
                         list(dayparts_df['end_hour'])
                         )
                  }

    months = list(seasons_dict)
    hours = list(range(1, 25))

    # APPLY TRANSFORMATION

    df_ts_template = pd.DataFrame(list(itertools.product(generation,
                                                         months,
                                                         hours,
                                                         years)
                                       ),
                                  columns=['TECHNOLOGY',
                                           'MONTH',
                                           'HOUR',
                                           'YEAR']
                                  )
    
    df_ts_template = df_ts_template.sort_values(by=['TECHNOLOGY', 'YEAR'])
    df_ts_template['SEASON'] = df_ts_template['MONTH'].map(seasons_dict)
    df_ts_template['DAYS'] = df_ts_template['SEASON'].map(days_dict)
    df_ts_template['YEAR'] = df_ts_template['YEAR'].astype(int)
    df_ts_template = powerplant_filter(df_ts_template)

    for daypart in dayparts_dict:
        if dayparts_dict[daypart][0] > dayparts_dict[daypart][1]: # loops over 24hrs
            df_ts_template.loc[(df_ts_template['HOUR'] >= dayparts_dict[daypart][0]) |
                          (df_ts_template['HOUR'] < dayparts_dict[daypart][1]),
                          'DAYPART'] = daypart
        else:
            df_ts_template.loc[(df_ts_template['HOUR'] >= dayparts_dict[daypart][0]) &
                      (df_ts_template['HOUR'] < dayparts_dict[daypart][1]),
                      'DAYPART'] = daypart

    df_ts_template = df_ts_template.drop_duplicates()
    df['SEASON'] = df['TIMESLICE'].str[0:2]
    df['DAYPART'] = df['TIMESLICE'].str[2:]
    df['YEAR'] = df['YEAR'].astype(int)
    df.drop(['REGION', 'FUEL', 'TIMESLICE'],
            axis=1,
            inplace=True)

    df = df.groupby(['LABEL', 'SEASON', 'DAYPART', 'YEAR'],
                    as_index=False)['VALUE'].sum()
    df = pd.merge(df,
                  df_ts_template,
                  how='left',
                  on=['LABEL', 'SEASON', 'DAYPART', 'YEAR']).dropna()
    df['HOUR_COUNT'] = df['DAYPART'].map(hours_dict)
    df['VALUE'] = (df['VALUE'].mul(1e6))/(df['DAYS']*df['HOUR_COUNT'].mul(3600))

    df = df.pivot_table(index=['MONTH', 'HOUR', 'YEAR'],
                        columns='LABEL',
                        values='VALUE',
                        aggfunc='mean').reset_index().fillna(0)
    df['MONTH'] = pd.Categorical(df['MONTH'],
                                 categories=months,
                                 ordered=True)
    df = df.sort_values(by=['MONTH', 'HOUR'])

    return df

def extract_codes_from_fuel(df: pd.DataFrame) -> pd.DataFrame:
    """Parses fuel codes into seperate columns
    
    Arguments:
        df: pd.DataFrame
            Dataframe with a FUEL column
    
    Returns:
        pd.DataFrame 
            Dataframe with FUEL code parsed into seperate columns
    """
    
    def sort_columns_fuel(df: pd.DataFrame) -> pd.DataFrame:
        for column in ("REGION", "COUNTRY", "REGION_CODE", "FUEL"):
            data = df.pop(column)
            df.insert(0, column, data)
        return df
    
    if "FUEL" not in df.columns:
        return df

    end_use = df.loc[df["FUEL"].str[-2:]=="02"]
    if not end_use.empty:
        end_use = end_use.drop(columns=["REGION"])
        end_use["F"] = end_use["FUEL"].str[0:3]
        end_use["REGION_CODE"] = end_use["FUEL"].str[3:8]
        end_use["COUNTRY"] = end_use["FUEL"].str[3:6]
        end_use["REGION"] = end_use["FUEL"].str[6:8]
        end_use = end_use.drop(columns=["FUEL"])
        end_use = end_use.rename(columns={"F":"FUEL"})
        end_use = sort_columns_fuel(end_use)
    else:
        end_use = end_use.drop(columns=["REGION", "TECHNOLOGY"])
        columns = list(df)
        columns.insert(0, "REGION")
        columns.insert(0, "COUNTRY")
        columns.insert(0, "REGION_CODE")
        columns.insert(0, "FUEL")
        end_use = pd.DataFrame(columns=columns)
        
    return end_use

def apply_timeshift(x, timeshift):
    """Applies timeshift to organize dayparts.
    
    Arguments:
        x = Value between 0-24
        timeshift = value offset from UTC (-11 -> +12)"""

    x += timeshift
    if x > 23:
        return x - 24
    elif x < 0:
        return x + 24
    else:
        return x