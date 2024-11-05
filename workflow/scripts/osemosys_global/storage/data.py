"""Functions to extract and format relevent data for storage."""
import pandas as pd

def get_years(start: int, end: int) -> range:
    return range(start, end + 1)

def get_years_build_rates(row):
    row['YEAR'] = get_years(row['START_YEAR'], row['END_YEAR'])
    
    return row

def format_gesdb_data(gesdb_data, 
                      gesdb_regional_mapping,
                      gesdb_tech_map,
                      inactive_vars,
                      active_vars,
                      new_vars
                      ):
    
    # Split 'storage_device' column to extract regional info
    data_storage_device = pd.json_normalize(gesdb_data['Storage Device'])
    
    data = pd.merge(gesdb_data, data_storage_device, 
                    left_index = True, right_index = True)

    data['Technology Mid-Type'] = data['Technology Mid-Type'
                                       ].fillna('null').replace('', 'null')
    data['Status'] = data['Status'].fillna('null')
    
    # Get unique region combinations to allow for OG mapping
    data['region'] = data['Country'].fillna('') + ';' + data[
        'State/Province/Territory'].fillna('')
    
    gesdb_regional_mapping['region'] = gesdb_regional_mapping[
        'Country'].fillna('') + ';' + gesdb_regional_mapping[
            'State/Province/Territory'].fillna('')
    
    # Check if all regions are defined in mapping file
    for region in data['region'].unique():
        if gesdb_regional_mapping['region'].str.contains(region).any() == False:
            raise Exception(
                f'{region} in storage dataset but not in GESDB_region_mapping.csv')
            
    gesdb_regional_mapping = gesdb_regional_mapping[gesdb_regional_mapping['node'].notna()]
    gesdb_regional_mapping = dict(zip(gesdb_regional_mapping['region'], gesdb_regional_mapping['node']))

    #Only keep entries with OG region
    data = data[data['region'].isin(list(gesdb_regional_mapping.keys()))]
    
    for tech in data['Technology Mid-Type'].unique():
        if tech not in list(gesdb_tech_map.keys()):
            raise Exception(f'{tech} technology in storage dataset but not in GESDB_TECH_MAP dict')
            
    for var in data['Status'].unique():
        if var not in inactive_vars and var not in active_vars and var not in new_vars:
            raise Exception(f'{var} status in storage dataset but not in ACTIVE_VARS or INACTIVE_VARS')
            
    # Add technology groups
    data['tech'] = data['Technology Mid-Type'].map(gesdb_tech_map)
    
    # Filter out entries without defined power rating
    data = data[data['Rated Power (kW)'] > 0]
    
    # Remove plants that are retired
    data = data.loc[~data['Status'].isin(inactive_vars)]
    
    # Add nodes
    data['REGION'] = data['region'].map(gesdb_regional_mapping)
    
    # Add OSeMOSYS Global naming convention
    data['STORAGE'] = data['tech'] + data['REGION'] + '01'

    usecols = ['Status', 'Rated Power (kW)', 'Storage Capacity (kWh)_x', 
               'REGION', 'Commissioned Date', 'tech', 'STORAGE']

    # Filter and rename relevant columns
    data = data[usecols].rename(columns = {'Rated Power (kW)' : 'rated_power', 
                                            'Storage Capacity (kWh)_x' : 'storage_capacity'})         
    
    return data, gesdb_regional_mapping