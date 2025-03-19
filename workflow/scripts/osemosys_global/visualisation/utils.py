"""Module for utility plotting functions."""

import pandas as pd
from typing import Dict, List, Union, Tuple
import itertools
from pathlib import Path
from constants import DAYS_PER_MONTH, MONTH_NAMES
from osemosys_global.utils import apply_timeshift
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

import logging 
logger = logging.getLogger(__name__)

def get_years(start: int, end: int) -> range:
    return range(start, end + 1)

def read_csv(dirpath: str) -> Dict[str,pd.DataFrame]:
    """Reads in CSVs folder
    
    Replace with ReadCSV.read() from otoole v1.0
    """
    data = {}
    files = [Path(x) for x in Path(dirpath).iterdir()]
    for f in files:
        data[f.stem] = pd.read_csv(f)
    return data

def filter_transmission_techs(df: pd.DataFrame, column_name: str = "TECHNOLOGY") -> pd.DataFrame:
    """Filters out only transmission technologies
    
    Arguments: 
        df: pd.DataFrame
            otoole formatted dataframe 
        column_name: str
            Column name to filter on 
            
    Returns: 
        pd.DataFrame
    """
    return df.loc[df[column_name].str.startswith("TRN")].reset_index(drop=True)

def get_color_codes(color_codes) -> Dict:
    """Get color naming dictionary.
    
    Return:
        Dictionary with tech and color name map
    """

    # Get colour mapping dictionary
    color_dict = dict([(n, c) for n, c
                   in zip(color_codes.tech_id,
                          color_codes.colour)])
    return color_dict

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
        filtered_df['LABEL'] = filtered_df['TYPE']
    else:
        filtered_df['LABEL'] = filtered_df['TYPE']
    
    filtered_df.drop(['TECHNOLOGY', 'TYPE', 'COUNTRY'],
            axis=1,
            inplace=True)
    return filtered_df

def transform_ts(seasons, 
                 dayparts,
                 timeshift,
                 start_year,
                 end_year,
                 data:Dict[str, pd.DataFrame], 
                 df:pd.DataFrame) -> pd.DataFrame:
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

    seasonsData = []

    for s, months in seasons.items():
        for month in months:
            seasonsData.append([month, s]) 
    seasons_df = pd.DataFrame(seasonsData, 
                              columns=['month', 'season'])
    seasons_df = seasons_df.sort_values(by=['month']).reset_index(drop=True)

    daypartData = []
    for dp, hr in dayparts.items():
        daypartData.append([dp, hr[0], hr[1]])
    dayparts_df = pd.DataFrame(daypartData,
                               columns=['daypart', 'start_hour', 'end_hour'])

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

    years = get_years(start_year, end_year[0])

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

def get_map(extent:List[Union[int,float]] = None) -> Tuple[plt.figure, plt.axes]: 
    """Sets up map to plot on
    
    Arguments: 
        extent:List[Union[int,float]]
            Lats and lons to plot formatted as 
            [lon_min, lon_max, lat_min, lat_max]
    """
    
    mrc = ccrs.Mercator()
    fig, ax = plt.subplots(subplot_kw={"projection": mrc})
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='lightgrey'))
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='black', linewidth = 0.3, facecolor='#46bcec'))
    ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth = 0.3)
    
    if extent:
        ax.set_extent(extent)
    
    return fig, ax

def plot_map_trn_line(ax: plt.axes, x1:float, y1:float, x2:float, y2:float, width:float) -> plt.axes:
    """Adds a transmission line to a map axes"""
    
    return ax.plot(
        [x1, x2], 
        [y1, y2],
        linewidth = width,
        color = 'crimson',
        transform=ccrs.PlateCarree(),
    )

def plot_map_text(ax: plt.axes, x:float, y:float, text:Union[int,float,str]) -> plt.axes:
    """Adds text to the map plot"""
    
    return ax.text(
        x, 
        y, 
        text, 
        transform = ccrs.PlateCarree(), 
        fontsize = 6, 
        fontweight= 'normal', 
        ha = 'right', 
        va = 'top', 
        alpha = .8
    )