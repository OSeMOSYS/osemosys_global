from dash import Dash, dcc, html
import geopandas as gpd
import numpy as np
import plotly.express as px
from osemosys_global.dashboard.components import ids
import osemosys_global.dashboard.constants as const
import logging 
logger = logging.getLogger(__name__)

def map_theme(**kwargs) -> html.Div:
    
    if "style" in kwargs:
        style = kwargs["style"]
    else:
        style = {}
    
    return html.Div(
        children=[
            html.H4("Map Theme"),
            dcc.Dropdown(
                id=ids.MAP_THEME_DROPDOWN,
                options=[
                    {"label":"Carto Darkmatter","value":"carto-darkmatter"},
                    {"label":"Carto Positron","value":"carto-positron"},
                    {"label":"Open Street Map","value":"open-street-map"},
                    {"label":"Stamen Terrain","value":"stamen-terrain"},
                    {"label":"Stamen Toner","value":"stamen-toner"},
                    {"label":"Stamen Watercolour","value":"stamen-watercolor"},
                    {"label":"White","value":"white-bg"}
                ],
                value=const._MAP_THEME,
                clearable=False,
                style={"vertical-align":"top"}
            )
        ],
        style=style,
    )

def transmission_lines(**kwargs) -> html.Div:
    
    if "style" in kwargs:
        style = kwargs["style"]
    else:
        style = {}
    
    return html.Div(
        children=[
            html.H4("Transmission Lines"),
            dcc.RadioItems(
                id=ids.NODE_TRANSMISSION_RADIO_BUTTON,
                options=["Show", "Hide"],
                value=const._TRANSMISSION_LINE,
                inline=True,
                style={"vertical-align":"top"}
            )
        ],
        style=style
    )
    
def node_locations(**kwargs) -> html.Div:
    
    if "style" in kwargs:
        style = kwargs["style"]
    else:
        style = {}
    
    return html.Div(
        children=[
            html.H4("Node Locations"),
            dcc.RadioItems(
                id=ids.NODE_LOCATION_RADIO_BUTTON,
                options=["Centroid", "Demand Center"],
                value=const._NODE_LOCATION,
                inline=True,
                style={"vertical-align":"top"}
            )
        ],
        style=style
    )

def map_size(**kwargs) -> html.Div:
    
    if "style" in kwargs:
        style = kwargs["style"]
    else:
        style = {}
    
    return html.Div(
        children=[
            dcc.Slider(
                id=ids.MAP_SIZE_SLIDER,
                min=1, 
                max=20, 
                step=1, 
                value=const._MAP_ELEMENTS_SIZE,
                included=False
            )
        ],
        style=style
    )

def plot_map(
    regions: list[str], 
    size: int, 
    show_lines: str, 
    map_theme: str, 
    nodes: gpd.GeoDataFrame,
    lines: gpd.GeoDataFrame
) -> html.Div:
    """Logic to plot map"""
    def include_connections(x:str, regions: list[str]) -> str:
        if (x[3:8] in regions) and (x[8:] in regions):
            return "Included"
        else:
            return "Not Included"
    
    if show_lines == "Hide":

        nodes["CATEGORY"] = nodes["NODE"].map(lambda x: "Included" if x in regions else "Not Included")
        nodes["SIZE"] = 1

        fig = px.scatter_mapbox(
            nodes,
            lat=nodes.geometry.y,
            lon=nodes.geometry.x,
            hover_name="NODE",
            hover_data={
                "SIZE":False,
                "CATEGORY":False
            },
            zoom=1,
            color=nodes["CATEGORY"],
            size=nodes["SIZE"],
            size_max=size,
            color_discrete_map={"Included":"#2e8b57", "Not Included":"#ff4500"},
            height=600
        )
    
    else:
        
        lines["CATEGORY"] = lines["TECHNOLOGY"].map(lambda x: include_connections(x, regions))
        lines["SIZE"] = 1
        
        lats_include = []
        lons_include = []
        names_include = []
        sizes_include = []
        
        lats_exclude = []
        lons_exclude = []
        names_exclude = []
        sizes_exclude = []

        for feature, name in zip(lines["geometry"], lines["TECHNOLOGY"]):
            
            x,y = feature.xy
            category = include_connections(name, regions)
            if category == "Included":
                lats_include = np.append(lats_include, y)
                lons_include = np.append(lons_include, x)
                names_include = np.append(names_include, [name]*len(y))
                sizes_include = np.append(sizes_include, [1]*len(y))

                # to create empty lines when switching between connectors
                lats_include = np.append(lats_include, None)
                lons_include = np.append(lons_include, None)
                names_include = np.append(names_include, None)
                sizes_include = np.append(sizes_include, None)
            else:
                lats_exclude = np.append(lats_exclude, y)
                lons_exclude = np.append(lons_exclude, x)
                names_exclude = np.append(names_exclude, [name]*len(y))
                sizes_exclude = np.append(sizes_exclude, [1]*len(y))

                # to create empty lines when switching between connectors
                lats_exclude = np.append(lats_exclude, None)
                lons_exclude = np.append(lons_exclude, None)
                names_exclude = np.append(names_exclude, None)
                sizes_exclude = np.append(sizes_exclude, None)

        fig = px.line_mapbox(
            lat=lats_include,
            lon=lons_include,
            hover_name=names_include,
            zoom=1,
            height=600
        )
        fig.add_trace(
            px.line_mapbox(
                lat=lats_exclude,
                lon=lons_exclude,
                hover_name=names_exclude,
                zoom=1,
                height=1000
            ).data[0]
        )
        fig.update_traces(line=dict(width=(size/4), color="#2e8b57"), selector=0)
        fig.update_traces(line=dict(width=(size/4), color="#ff4500"), selector=1)

    fig.update_layout(mapbox_style=map_theme)
    return html.Div(children=[dcc.Graph(figure=fig)], id=ids.MAP)
