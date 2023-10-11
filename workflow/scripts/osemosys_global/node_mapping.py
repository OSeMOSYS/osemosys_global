import pandas as pd
import os
import sys
import yaml
import geopandas as gpd
import matplotlib.pyplot as plt

from OPG_configuration import ConfigFile, ConfigPaths
pd.set_option('mode.chained_assignment', None)


def main():
    
    # Config paths
    config_paths = ConfigPaths()
    input_data_dir = config_paths.input_data_dir

    # Mapping admin_1 to custom nodes
    asean_custom_nodes = pd.read_csv(os.path.join(input_data_dir,
                                                  'asean_nodes_mapping_24.csv'))
    print(asean_custom_nodes)
    
    feo_global_nodes = pd.read_csv(os.path.join(input_data_dir,
                                               'feo_global_node_mapping_final.csv'))
    
    # Read geopackage with admin 1 at global coverage 
    gdf = gpd.read_file(os.path.join(input_data_dir,
                                     'gadm_410_admin_1_with_node_codes.gpkg'))  
    
    # Filter to show only the 10 ASEAN countries
    asean_list = ('KHM',
                  'LAO',
                  'VNM',
                  'MMR',
                  'THA',
                  'MYS',
                  'SGP',
                  'BRN',
                  'PHL',
                  'IDN')
    
    feo_global_list = tuple(feo_global_nodes['admin_0_code'].unique())
    print(feo_global_list)
    
    gdf = gdf.loc[gdf['node_code'].str.startswith(feo_global_list)]
    gdf = pd.merge(gdf, feo_global_nodes,
                   how='left',
                   on=['node_code'])
    
    print(gdf.head())
    gdf_custom = gdf.dissolve(by='custom_node_code')
    gdf_custom.reset_index(inplace=True)
    ax = gdf_custom.plot(column=gdf_custom['custom_node_code'].str[:3],
                         cmap='tab20')
    ax.set_axis_off()
    gdf_custom.apply(lambda x: ax.annotate(text=x['custom_node_code'],
                                           xy=x.geometry.centroid.coords[0],
                                           ha='center',
                                           fontsize=8,
                                           font='Arial'),
                     axis=1)
    plt.figure(figsize=(6,10))
    plt.show()
    
    gdf_custom.to_file(os.path.join(input_data_dir,
                                    'feo_global_custom_nodes.gpkg'),
                       driver='GPKG')
    
if __name__ == '__main__':
    main()
