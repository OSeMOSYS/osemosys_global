"""Function to calculate residual capacity for transmission technologies."""
import pandas as pd
import numpy as np

from data import(
    get_years, 
    format_gesdb_data
    )

def res_capacity_storage(gesdb_data, 
                         gesdb_regional_mapping,
                         res_cap_base,
                         op_life_dict,
                         storage_param,
                         gesdb_tech_map,
                         duration_type,
                         build_year,
                         retirement_year,
                         inactive_vars,
                         active_vars,
                         new_vars,
                         start_year,
                         end_year,
                         region_name
                         ):
    
    data, regional_mapping = format_gesdb_data(
        gesdb_data, 
        gesdb_regional_mapping,
        gesdb_tech_map,
        inactive_vars,
        active_vars,
        new_vars
        )
    
    duration_dict = {}
    
    years = get_years(start_year, end_year)
    
    # Set baseline technology capital costs and duration as defined in the config file.
    for tech, tech_params in storage_param.items():
        duration_dict[tech] = tech_params[4]
    
    # Set historic duration if 'Historic' is chosen for 'duration_type'.
    if duration_type == 'Historical':
        # Set total hours for entries that have both variables defined.
        data['hours'] = data['storage_capacity'] / data['rated_power']
        
        # Calculte average storage hours per aggregrate technology group.
        avg_hours = data.groupby('tech')['hours'].mean().reset_index()
        avg_hours = dict(zip(avg_hours['tech'], avg_hours['hours']))
        
        # Assign average storage hours per aggregrrate technology group if entry specific hours is not defined.
        data['hours'] = np.where((data['hours'] == 0 )| (np.isnan(data['hours'])), data['tech'].map(duration_dict), data['hours'])
       
    # Set user defined storage duration if 'Default' is chosen for 'duration_type'.
    else:
        data['hours'] = data['tech'].map(duration_dict)

    data['storage_capacity'] = data['rated_power'] * data['hours']

    # Calculate commissioning year for operational plants without commissioning year entry
    data['Commissioned Date'] = data['Commissioned Date'].str[-4:].astype('Int32')

    avg_date = round(data.groupby('tech')['Commissioned Date'].mean().reset_index(),0)
    avg_date = dict(zip(avg_date['tech'], avg_date['Commissioned Date']))

    data['Commissioned Date'] = np.where(data['Commissioned Date'].isna() & data['Status'].isin(active_vars), 
                                         data['tech'].map(avg_date), data['Commissioned Date'])

    # Fill commissioning year for new plants without commissioning year entry
    data['Commissioned Date'] = np.where(data['Commissioned Date'].isna() & data['Status'].isin(new_vars), 
                                         build_year, data['Commissioned Date'])

    # Add retirement year for operational/new plants
    data['Retirement Date'] = np.where(data['Commissioned Date'] + data['tech'].map(op_life_dict) < start_year + 1, 
                                       retirement_year, data['Commissioned Date'] + data['tech'].map(op_life_dict))

    # Aggregate entries for existing storage plants.
    data_existing = data.loc[data['Commissioned Date'] <= start_year].groupby(['STORAGE']
                                ).sum(['rated_power', 'storage_capacity']).drop(columns = {'hours'}).reset_index()

    for tech in data['STORAGE'].unique():
        if data_existing.loc[data_existing['STORAGE'] == tech].empty:
            data_existing.loc[len(data_existing.index)] = [tech, 0, 0, '-', '-']

    data_existing['YEAR'] = start_year

    # Aggregate entries for new storage plants.
    data_new = data.loc[data['Commissioned Date'] > start_year]

    data_new = data.loc[data['Commissioned Date'] > start_year].groupby(['STORAGE', 'Commissioned Date', 'Retirement Date']
                                ).sum(['rated_power', 'storage_capacity']).drop(columns = {'hours'}).reset_index(
                                    ).rename(columns = {'Commissioned Date' : 'YEAR'})

    data_new['Status'] = 'new'

    # Aggregate entries for storage plants that retire before the end of the model horizon.
    data_retirement = data.loc[data['Retirement Date'] <= end_year].groupby(['STORAGE', 'Commissioned Date', 'Retirement Date']
                                ).sum(['rated_power', 'storage_capacity']).drop(columns = {'hours'}).reset_index(
                                    ).rename(columns = {'Retirement Date' : 'YEAR'})

    data_retirement['Status'] = 'retired'

    # Combine dfs new and retirement dfs.
    data_change = pd.concat([data_new, data_retirement]
                            ).sort_values(['YEAR']).reset_index(drop = True
                                                                ).drop(columns = {'Commissioned Date', 'Retirement Date'})

    # Calculate residual capacity in terms of storage capacity (PJ) and discharge rates (PJ/hr)
    residual = data_existing.copy()

    # Add new capacities/retired capacities by horizon year
    for entry in data_change.index:
        change = data_change.loc[data_change.index == entry].reset_index(drop = True)
        # Pull the old/existing data from the storage entry that requires change.
        old = residual.loc[residual['STORAGE'] == change['STORAGE'].iloc[0]].iloc[[-1]].reset_index(drop = True)
        
        # If the required change is a new entry sum the entry data and concat to residuals.
        if change['Status'].iloc[0] == 'new':
            change['rated_power'] = change['rated_power'] + old['rated_power']
            change['storage_capacity'] = change['storage_capacity'] + old['storage_capacity']
            residual = pd.concat([residual, change], join = 'inner')

        # If the required change is a new entry substract the entry data and concat to residuals.
        if change['Status'].iloc[0] == 'retired':
            change['rated_power'] = old['rated_power'] - change['rated_power']
            change['storage_capacity'] = old['storage_capacity'] - change['storage_capacity']
            residual = pd.concat([residual, change], join = 'inner')
            
    for tech in residual['STORAGE'].unique():
        for year in years:
            if residual.loc[(residual['STORAGE'] == tech) & 
                            (residual['YEAR'] == year)].empty:
                
                new_row = residual.loc[(residual['STORAGE'] == tech) & 
                                       (residual['YEAR'] == (year - 1))]
                new_row.loc[new_row['YEAR'] == (year - 1), 'YEAR'] = year
                
                residual = pd.concat([residual, new_row])           

    residual['REGION'] = region_name
    residual = residual.drop_duplicates(subset=['REGION', 'STORAGE', 'YEAR']
                                        , keep = 'last').set_index(['REGION', 'STORAGE', 'YEAR'])
    
    # Pull power rating data for ResidualCapacity.csv
    residual_capacity = residual[['rated_power']].reset_index(
        drop = False).rename(columns = {'rated_power' : 'VALUE', 'STORAGE' : 'TECHNOLOGY'})
    
    # Convert from kW to GW
    residual_capacity['VALUE'] = round(residual_capacity['VALUE'] / 1000000, 4)
    residual_capacity['TECHNOLOGY'] = 'PWR' + residual_capacity['TECHNOLOGY']
    
    # Pull storage capacity data for ResidualStorageCapacity.csv
    residual_storage_capacity = residual[['storage_capacity']].reset_index(
        drop = False).rename(columns = {'storage_capacity' : 'VALUE'})
    
    # Convert from kWh to PJ
    residual_storage_capacity['VALUE'] = round(residual_storage_capacity['VALUE'
                                                                         ] / 277777777.77778, 5)
 
    # Combine residual capacity df's.
    residual_capacity = pd.concat([res_cap_base, residual_capacity]).reset_index(drop = True)

    return residual_capacity, residual_storage_capacity