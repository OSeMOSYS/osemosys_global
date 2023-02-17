from pathlib import Path
from dash import Dash
import i18n 
from ..components.layout import create_layout
from data.loader import load_osemosys_data, load_node_data, load_line_data
from ..utils import load_connections, load_connections_2, seperate_connectors
from ..OPG_configuration import ConfigFile, ConfigPaths

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

LOCALE = "en"

def main() -> None:
    
    # set mapping options
    i18n.set("locale",LOCALE)
    i18n.load_path.append("./locale")
    
    # read in data
    config = ConfigPaths()
    osemosys_data = load_osemosys_data(config.output_data_dir)
    node_data = load_node_data("./nodes.csv")
    
    # set up node connection data 
    connections = load_line_data("./lines.csv")
    # connections = load_connections(connections, node_data)
    # connections = seperate_connectors(connections)
    connections = load_connections_2(connections, node_data)
    
    # create app
    app = Dash(external_stylesheets=external_stylesheets)
    app.title = i18n.t("general.app_title")
    app.layout = create_layout(
        app, 
        osemosys_data, 
        node_data, 
        connections)
    app.run()

if __name__ == "__main__":
    main()

