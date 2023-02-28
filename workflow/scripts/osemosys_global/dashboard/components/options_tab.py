"""Module to hold components of options tab"""

from typing import List
import pandas as pd
from dash import Dash, dcc, html
from typing import Dict
from osemosys_global.dashboard.components import ids
import osemosys_global.dashboard.constants as const
from osemosys_global.dashboard.utils import to_dropdown_options

def plot_theme() -> html.Div:
    return html.Div(
        children=[
            html.H3("Plot Theme"),
            dcc.Dropdown(
                id=ids.PLOT_THEME_DROPDOWN,
                options=[
                    {"label":"Plotly Standard","value":"plotly"}, 
                    {"label":"Plotly White","value":"plotly_white"}, 
                    {"label":"Plotly Dark","value":"plotly_dark"}, 
                    {"label":"ggplot","value":"ggplot2"}, 
                    {"label":"Seaborn","value":"seaborn"}, 
                    {"label":"Simple White","value":"simple_white"}
                ],
                value=const._PLOT_THEME,
                clearable=False,
            ),
        ]
    )
    
def plotting_level() -> html.Div:
    return html.Div(
        children=[
            html.H3("Plotting Level"),
            dcc.RadioItems(
                id=ids.REGION_COUNTRY_RADIO_BUTTON,
                options=["System", "Country", "Region"],
                value= const._GEOGRAPHIC_SCOPE,
                inline=True
            ),
        ]
    )
    
def country_dropdown(countries: List[str]) -> html.Div:
    return html.Div(
        children=[
            html.H3("Countries to Include"),
            dcc.Dropdown(
                id=ids.COUNTRY_DROPDOWN,
                options=to_dropdown_options(countries),
                value=countries,
                multi=True,
                persistence=True
            ),
            html.Button(
                children=["Select All"],
                id=ids.SELECT_ALL_COUNTRIES_BUTTON,
                n_clicks=0
            ),   
        ]
    )
    
def region_dropdown(regions: List[str]) -> html.Div:
    return html.Div(
        children=[
            html.H3("Regions to Include"),
            dcc.Dropdown(
                id=ids.REGION_DROPDOWN,
                options=to_dropdown_options(regions),
                value=regions,
                multi=True,
                persistence=True
            ),
            html.Button(
                children=["Select All"],
                id=ids.SELECT_ALL_REGIONS_BUTTON,
                n_clicks=0
            )
        ]
    )
