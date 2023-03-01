"""Holds components for result data tab"""

from dash import Dash, html, dcc
import osemosys_global.dashboard.constants as const
import osemosys_global.dashboard.components.ids as ids
import pandas as pd
from typing import Dict, List
from osemosys_global.dashboard.utils import (
    get_unique_techs, 
    get_unique_fuels,
    create_dropdown_options
)

def variable_dropdown(**kwargs) -> html.Div:
    """Selects model variable"""
    
    if "style" in kwargs:
        style = kwargs["style"]
    else:
        style = {}
    
    return html.Div(
        children=[
            html.H3("Variable"),
            dcc.Dropdown(
                id=ids.RESULT_DATA_DROPDOWN,
                options = [{"label":const.RESULT_CONFIG[x]["nicename"], "value":x} for x in const.RESULT_CONFIG],
                value=const._RESULT_DATA_CHOICE,
                multi=False,
                clearable=False
            )
        ],
        style=style
    )

def year_selector(years: List[int]) -> html.Div:
    return html.Div(
        children=[
            dcc.Slider(
                id=ids.RESULT_YEAR_SELECTOR,
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
                id=ids.RESULT_YEAR_SLIDER,
                min=years[0], 
                max=years[-1], 
                step=10, 
                value=[years[0], years[-1]],
                marks={i:f"{int(i)}" for i in years}
            )
        ]
    )

def tech_dropdown(input_data: Dict[str,pd.DataFrame], parameter: str, geographic_scope: str, **kwargs) -> html.Div:
    """Dropdown for the tech/fuel selection
    
    Arguments:
        input_data: Dict[str,pd.DataFrame]
            Internal datastore
        parameter: str
            Parameter to write options for 
        geographic_scope: str
            Plotting level 
    """
    if const.RESULT_CONFIG[parameter]["groupby"] == "TECHNOLOGY":
        options = get_unique_techs(input_data[parameter])
    else:
        options = get_unique_fuels(input_data[parameter])
        
    label_values = create_dropdown_options(options)
    
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
                id=ids.RESULT_DATA_TECH_DROPDOWN,
                options = label_values,
                value=label_values[0]["value"],
                multi=False,
                clearable=False
            )
        ],
        style=style
    )
    
