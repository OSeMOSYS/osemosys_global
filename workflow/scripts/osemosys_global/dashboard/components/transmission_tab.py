"""Module to hold transmission result components"""

from dash import html, dcc
import osemosys_global.dashboard.components.ids as ids
import osemosys_global.dashboard.constants as const
from osemosys_global.dashboard.utils import get_transmission_techs, plot_by_system
from typing import List, Dict
import pandas as pd

def line_dropdown(lines: List[str], **kwargs) -> html.Div:
    """Selects model variable"""
    
    if "style" in kwargs:
        style = kwargs["style"]
    else:
        style = {}
    
    options = [{"label":x, "value":x} for x in lines]
    options.insert(0,{"label":"All", "value":"all"})
    
    return html.Div(
        children=[
            html.H3("Transmission Line"),
            dcc.Dropdown(
                id=ids.TRANSMISSION_LINE_DROPDOWN,
                options = options,
                value=lines[0],
                multi=False,
                clearable=False
            )
        ],
        style=style
    )
    
def year_selector(years: List[int]) -> html.Div:
    """Single year selector"""
    return html.Div(
        children=[
            dcc.Slider(
                id=ids.TRANSMISSION_YEAR_SELECTOR,
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
    """Multi-year selector"""
    return html.Div(
        children=[
            dcc.RangeSlider(
                id=ids.TRANSMISSION_YEAR_SLIDER,
                min=years[0], 
                max=years[-1], 
                step=10, 
                value=[years[0], years[-1]],
                marks={i:f"{int(i)}" for i in years}
            )
        ]
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
                id=ids.TRANSMISSION_DATA_DROPDOWN,
                options = [{"label":const.TRANSMISSION_CONFIG[x]["nicename"], "value":x} for x in const.TRANSMISSION_CONFIG],
                value=const._TRANSMISSION_DATA_CHOICE,
                multi=False,
                clearable=False
            )
        ],
        style=style
    )
    
def plot_transmission_data(
    data: pd.DataFrame,
    line: str,
    years: List[int], 
    parameter: str,
    plot_theme: str,
    plot_type: str,
    config: Dict[str,Dict],
    div_id: str,
) -> html.Div:
    """Generic plotting function for transmission data"""
    
    # parse input data codes 
    data = get_transmission_techs(data)
    if line != "all":
        data = data.loc[data["TECHNOLOGY"] == line].reset_index(drop=True)
        
    # temporal filtering 
    if len(years) == 1:
        filtered = data.loc[data["YEAR"] == years[0]]
    else:
        filtered = data.loc[data["YEAR"].isin(range(years[0],(years[-1]+1)))]

    if filtered.shape[0] == 0:
        return html.Div("No Data Selected")
    
    # aggregate data
    groupby_column = config[parameter]["xaxis"]
    grouped = filtered.groupby(by=[groupby_column]).sum(numeric_only=True).reset_index()
    fig = plot_by_system(
        df=grouped,
        plot_type=plot_type,
        x=config[parameter]["xaxis"],
        y="VALUE",
        plot_theme=plot_theme,
        labels={"VALUE":config[parameter]["ylabel"]}
    )

    return html.Div(children=[dcc.Graph(figure=fig)], id=div_id)