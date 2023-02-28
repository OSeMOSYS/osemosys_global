"""Utility functions for dashboard"""

import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString
from osemosys_global.visualisation.utils import load_node_data_demand_center, load_node_data_centroid, load_line_data
import plotly.express as px
from typing import Union, List

def add_pts_to_line(x1: float, y1: float, x2: float, y2: float, n_points: int, name: str) -> list[str, LineString]:
    """Generate a straight line of coordinates between two input 2D coordinates"""
    
    x_coords = np.linspace(x1, x2, n_points)
    y_coords = np.linspace(y1, y2, n_points)
    return [name, LineString(zip(x_coords, y_coords))]

def geolocate_nodes(node_data_file: str, centroid: bool) -> gpd.GeoDataFrame:
    """Geolocates node data
    
    Arguments:
        node_data_file: str
            Path to datafile that holds node data  
        centroid: bool
            Get centroid data, else, demand center data
            
    Returns: 
        GeoDataframe with nodes and lats and lons
    """

    if centroid:
        nodes = load_node_data_centroid(node_data_file)
    else:
        nodes = load_node_data_demand_center(node_data_file)
        
    return gpd.GeoDataFrame(
        nodes, geometry=gpd.points_from_xy(nodes["LONGITUDE"], nodes["LATITUDE"])
    )
    
def geolocate_lines(cost_line_expansion_xlsx: str, nodes: gpd.GeoDataFrame, n_points: int = 20):
    """Gets transmission lines with lats and lons
    
    Arguments:
        cost_line_expansion_xlsx: str
            Path to 'Costs Line expansion.xlsx' file  
        geolocated_nodes: gpd.GeoDataFrame
            Geolocated nodes, formatted from geolocate_nodes()
        n_points: int
            Number of points to represent on the output line
            
    Returns: 
        GeoDataframe with start/end nodes and lats/lons
    """
    
    trn = load_line_data(cost_line_expansion_xlsx)
    
    lats_lookup = dict(zip(nodes['NODE'], nodes['LATITUDE']))
    lons_lookup = dict(zip(nodes['NODE'], nodes['LONGITUDE']))
    data_to_expand = zip(trn["FROM"].to_list(), trn["TO"].to_list())
    
    data = []
    for from_region, to_region in data_to_expand:
        try:
            data.append(
                add_pts_to_line(
                    lons_lookup[from_region],
                    lats_lookup[from_region],
                    lons_lookup[to_region],
                    lats_lookup[to_region],
                    n_points,
                    f"TRN{from_region}{to_region}"
                )
            )
        except KeyError:
            print(f"Can not find {from_region} and/or {to_region} in node data")
    return gpd.GeoDataFrame(data, columns=["TECHNOLOGY", "geometry"])

def get_regions(
    input_data: dict[str, pd.DataFrame],
    countries_only: bool = False,
) -> list[str]:
    """Gets all region codes"""
    
    def filter_region_codes(df:pd.DataFrame) -> list:
        techs = df.loc[
            (df["VALUE"].str.startswith("PWR")) &
            ~(df["VALUE"].str.startswith("PWRTRN"))]
        regions = techs["VALUE"].str[-7:-2]
        return regions.drop_duplicates().to_list()
    
    region_data = filter_region_codes(input_data["TECHNOLOGY"])
    if countries_only:
        return list(set([x[0:3] for x in region_data]))
    else:
        return region_data
    
def to_dropdown_options(values: list[Union[str,int]]) -> list[dict[str,str]]:
    return [{"label": value, "value": value} for value in values]

def get_generation_techs(df: pd.DataFrame, column:str = "TECHNOLOGY") -> pd.DataFrame:
    """Filters out only power generation technologies"""
    return df.loc[
        (df[column].str.startswith("PWR")) & 
        ~(df[column].str.startswith("PWRTRN"))].reset_index(drop=True)
    
def parse_pwr_codes(df: pd.DataFrame) -> pd.DataFrame:
    """Strips out country/region/category tags from pwr codes"""

    def sort_columns(df: pd.DataFrame) -> pd.DataFrame:
        for column in ("REGION", "COUNTRY", "REGION_CODE", "CATEGORY"):
            data = df.pop(column)
            df.insert(0, column, data)
        return df
    
    if not df.empty:
        df = df.drop(columns=["REGION"])
        df["CATEGORY"] = df["TECHNOLOGY"].str[3:6]
        df["REGION_CODE"] = df["TECHNOLOGY"].str[6:11]
        df["COUNTRY"] = df["TECHNOLOGY"].str[6:9]
        df["REGION"] = df["TECHNOLOGY"].str[9:11]
        df = df.drop(columns=["TECHNOLOGY"])
        df = sort_columns(df)
    else:
        df = df.drop(columns=["REGION", "TECHNOLOGY"])
        columns = list(df)
        columns.insert(0, "CATEGORY")
        columns.insert(0, "REGION_CODE")
        columns.insert(0, "COUNTRY")
        columns.insert(0, "REGION")
        df = pd.DataFrame(columns=columns)
    return df

def parse_fuel_codes(df: pd.DataFrame) -> pd.DataFrame:
    """Parses fuel codes into seperate columns"""
    def sort_columns(df: pd.DataFrame) -> pd.DataFrame:
        for column in ("REGION", "COUNTRY", "REGION_CODE", "CATEGORY"):
            data = df.pop(column)
            df.insert(0, column, data)
        return df
    
    df = df.loc[df["FUEL"].str[-2:]=="02"]
    if not df.empty:
        df = df.drop(columns=["REGION"])
        df["CATEGORY"] = df["FUEL"].str[0:3]
        df["REGION_CODE"] = df["FUEL"].str[3:8]
        df["COUNTRY"] = df["FUEL"].str[3:6]
        df["REGION"] = df["FUEL"].str[6:8]
        df = df.drop(columns=["FUEL"])
        df = sort_columns(df)
    else:
        df = df.drop(columns=["REGION", "TECHNOLOGY"])
        columns = list(df)
        columns.insert(0, "REGION")
        columns.insert(0, "COUNTRY")
        columns.insert(0, "REGION_CODE")
        columns.insert(0, "CATEGORY")
        df = pd.DataFrame(columns=columns)
    return df
    
def plot_by_region(
    df: pd.DataFrame, 
    plot_type: str, 
    x: str,
    y: str,
    color: str,
    plot_theme: str = "plotly",
    **kwargs
) -> Union[px.bar, px.line, px.area]:
    """Plots a chart based on geography filtered data"""
    
    if plot_type == "Bar":
        return px.bar(df, x=x, y=y, color=color, barmode='group', template=plot_theme, **kwargs)
    elif plot_type == "Line":
        return px.line(df, x=x, y=y, color=color, template=plot_theme, **kwargs)
    elif plot_type == "Area":
        return px.area(df, x=x, y=y, color=color, template=plot_theme, **kwargs)

def plot_by_system(
    df: pd.DataFrame, 
    plot_type: str, 
    x: str,
    y: str,
    plot_theme: str = "plotly",
    **kwargs,
) -> Union[px.bar, px.line, px.area]:
    """Plots a chart based on system wide data"""

    if plot_type == "Bar":
        return px.bar(df, x=x, y=y, barmode='group', template=plot_theme, **kwargs)
    elif plot_type == "Line":
        return px.line(df, x=x, y=y, template=plot_theme, **kwargs)
    elif plot_type == "Area":
        return px.area(df, x=x, y=y, template=plot_theme, **kwargs)
    
def get_unique_techs(df: pd.DataFrame) -> List[str]:
    """Gets unique 3 letter tech codes
    
    Arguments:
        df: pd.DataFrame
            DataFrame to search through 
            
    Returns:
        List[str]
            Unique codes without region identifier. 
    """
    techs = get_generation_techs(df)
    techs["CATEGORY"] = techs["TECHNOLOGY"].str[3:6]
    return techs["CATEGORY"].unique()

def get_unique_fuels(df: pd.DataFrame) -> List[str]:
    """Gets unique 3 letter tech codes
    
    Arguments:
        df: pd.DataFrame
            DataFrame to search through 
            
    Returns:
        List[str]
            Unique codes without region identifier. 
    """
    fuels = df["FUEL"]
    return fuels.str[0:3].unique()
