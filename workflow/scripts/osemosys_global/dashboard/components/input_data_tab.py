"""Functions for input data tab"""

from typing import List
import pandas as pd
from dash import Dash, dcc, html
from typing import Dict
from osemosys_global.visualisation.utils import get_color_codes
from osemosys_global.dashboard.components import ids
import osemosys_global.dashboard.constants as const
from osemosys_global.dashboard.utils import get_generation_techs, parse_pwr_codes, parse_fuel_codes, plot_by_region, plot_by_system, get_unique_techs, get_unique_fuels

def year_selector(years: List[int]) -> html.Div:
    return html.Div(
        children=[
            dcc.Slider(
                id=ids.YEAR_SELECTOR_INPUT_DATA,
                min=years[0], 
                max=years[-1], 
                step=1, 
                value=years[0],
                marks={i:f"{int(i)}" for i in years},
                included=False
            )
        ]
    )
    
def year_slider(years: List[int]) -> html.Div:
    return html.Div(
        children=[
            dcc.RangeSlider(
                id=ids.YEAR_SLIDER_INPUT_DATA,
                min=years[0], 
                max=years[-1], 
                step=10, 
                value=[years[0], years[-1]],
                marks={i:f"{int(i)}" for i in years}
            )
        ]
    )

def input_data_dropdown(**kwargs) -> html.Div:
    
    if "style" in kwargs:
        style = kwargs["style"]
    else:
        style = {}
    
    return html.Div(
        children=[
            html.H3("Parameter"),
            dcc.Dropdown(
                id=ids.INPUT_DATA_DROPDOWN,
                options = [{"label":const.PARAM_CONFIG[x]["nicename"], "value":x} for x in const.PARAM_CONFIG],
                value=const._INPUT_DATA_CHOICE,
                multi=False,
                clearable=False
            )
        ],
        style=style
    )

def plot_type_dropdown(**kwargs) -> html.Div:
    
    if "style" in kwargs:
        style = kwargs["style"]
    else:
        style = {}
    
    return html.Div(
        children=[
            html.H3("Plot Type"),
            dcc.Dropdown(
                id=ids.PLOT_TYPE_DROPDOWN,
                options=["Bar", "Line", "Area"],
                value=const._PLOT_TYPE,
                clearable=False
            )
        ],
        style=style
    )

def tech_filter_dropdown(input_data: Dict[str,pd.DataFrame], parameter: str, geographic_scope: str, **kwargs) -> html.Div:
    """Dropdown for the tech/fuel selection
    
    Arguments:
        input_data: Dict[str,pd.DataFrame]
            Internal datastore
        parameter: str
            Parameter to write options for 
        geographic_scope: str
            Plotting level 
    """
    if const.PARAM_CONFIG[parameter]["groupby"] == "TECHNOLOGY":
        options = get_unique_techs(input_data[parameter])
    else:
        options = get_unique_fuels(input_data[parameter])
        
    label_values = assign_label_values(options)
    
    if geographic_scope == "System":
        label_values.insert(0, {"label":"All", "value":"all"})
    
    if "style" in kwargs:
        style = kwargs["style"]
    else:
        style = {}
    
    return html.Div(
        children=[
            html.H3("Technology/Fuel"),
            dcc.Dropdown(
                id=ids.INPUT_DATA_TECH_DROPDOWN,
                options = label_values,
                value=label_values[0]["value"],
                multi=False,
                clearable=False
            )
        ],
        style=style
    )

def assign_label_values(options: List[str]) -> List[Dict[str,str]]:
    """Creates label name map for dropdown"""
    try:
        return [{"label":const.TECHS_CONFIG[x]["nicename"], "value":x} for x in options]
    except KeyError:
        return[{"label":x, "value":x} for x in options]

def filter_data(data: pd.DataFrame, region_filter_column: str, regions: List[str], years: List[int]) -> pd.DataFrame:
    """Filter data by geography and time

    Arguments:
        data: pd.DataFrame
            Data to filter
        region_filter_column: str, 
            Column to filter regions on. Must be in set ("COUNTRY", "REGION")
        regions: List[str]
            Regions to filter over
        years: List[int]
            years to filter over
    """
    
    if len(years) == 1:
        return data.loc[
                (data[region_filter_column].isin(regions)) &
                (data["YEAR"] == years[0])
            ]
    else:
        return data.loc[
                (data[region_filter_column].isin(regions)) &
                (data["YEAR"].isin(range(years[0],(years[-1]+1))))
            ]

def group_data(data: pd.DataFrame, groupby_columns: List[str], groupby_method: str) -> pd.DataFrame:
    """Aggregates data based on spatial scope
    
    Arguments:
        data: pd.DataFrame
            Data to filter
        groupby_columns = List[str]
            columns to group on 
    """
    if groupby_method == "sum":
        return data.groupby(by=groupby_columns).sum(numeric_only=True).reset_index()
    elif groupby_method == "mean":
        return data.groupby(by=groupby_columns).mean(numeric_only=True).reset_index()
    else:
        return pd.DataFrame()

def plot_input_data(
    input_data: Dict[str, pd.DataFrame],
    countries: list[str],
    regions: list[str],
    plot_theme: str,
    geographic_scope: str,
    years: list[int], 
    plot_type: str,
    parameter: str,
    tech_fuel_filter: str
) -> html.Div:
    """Generic function for plotting input data"""
    
    # parse input data codes 
    data = input_data[parameter]
    if const.PARAM_CONFIG[parameter]["groupby"] == "TECHNOLOGY":
        data = get_generation_techs(data)
        data = parse_pwr_codes(data)
        if tech_fuel_filter != "all":
            data = data.loc[data["CATEGORY"] == tech_fuel_filter]
    else:
        data = parse_fuel_codes(data)
        if tech_fuel_filter != "all":
            data = data.loc[data["CATEGORY"] == tech_fuel_filter]

    # determine filetering options 
    if geographic_scope == "Country":
        region_filter_column = "COUNTRY"
        region_filters = countries
        plot_color = "COUNTRY"
        groupby_columns = ["CATEGORY", region_filter_column, const.PARAM_CONFIG[parameter]["xaxis"]]
    elif geographic_scope == "Region":
        region_filter_column = "REGION_CODE"
        region_filters = regions
        plot_color = "REGION_CODE"
        groupby_columns = ["CATEGORY", region_filter_column, const.PARAM_CONFIG[parameter]["xaxis"]]
    else: # system level
        region_filter_column = "REGION_CODE"
        region_filters = regions
        if tech_fuel_filter != "all":
            groupby_columns = ["CATEGORY", const.PARAM_CONFIG[parameter]["xaxis"]]
            plot_color = const.PARAM_CONFIG[parameter]["groupby"]
            color_map = get_color_codes()
        else: 
            groupby_columns = ["CATEGORY", const.PARAM_CONFIG[parameter]["xaxis"]]
            plot_color = "CATEGORY"
            color_map = get_color_codes()
        
    groupby_method = const.PARAM_CONFIG[parameter]["groupby_method"]
        
    # get regional and temporal filtered dataframe 
    filtered_data = filter_data(
        data = data,
        region_filter_column = region_filter_column,
        regions = region_filters, 
        years = years
    )
    if filtered_data.shape[0] == 0:
        return html.Div("No Data Selected")
    
    # aggregate data
    grouped_data = group_data(data=filtered_data, groupby_columns=groupby_columns, groupby_method=groupby_method)
    
    # plot data
    if geographic_scope in ["Country", "Region"]:
        fig = plot_by_region(
            df=grouped_data,
            plot_type=plot_type,
            x=const.PARAM_CONFIG[parameter]["xaxis"],
            y="VALUE",
            color=plot_color,
            plot_theme=plot_theme,
            labels={"VALUE":const.PARAM_CONFIG[parameter]["ylabel"]}
        )
    else:
        if tech_fuel_filter != "all":
            fig = plot_by_system(
                df=grouped_data,
                plot_type=plot_type,
                x=const.PARAM_CONFIG[parameter]["xaxis"],
                y="VALUE",
                plot_theme=plot_theme,
                labels={"VALUE":const.PARAM_CONFIG[parameter]["ylabel"]}
            )
        else:
            fig = plot_by_system(
                df=grouped_data,
                plot_type=plot_type,
                x=const.PARAM_CONFIG[parameter]["xaxis"],
                y="VALUE",
                color=plot_color,
                color_discrete_map=color_map,
                plot_theme=plot_theme,
                labels={"VALUE":const.PARAM_CONFIG[parameter]["ylabel"]}
            )

    return html.Div(children=[dcc.Graph(figure=fig)], id=ids.INPUT_DATA_CHART)
