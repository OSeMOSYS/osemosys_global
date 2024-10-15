"""Functions to extract and format relevent data for tranmission."""
import pandas as pd
import geopy.distance

from sets import get_unique_technologies

def get_years(start: int, end: int) -> range:
    return range(start, end + 1)

def format_gtd_existing(df):
    cols = {'pathway' : 'TECHNOLOGY',
            'from_region' : 'From',
            'to_region' : 'To', 
            'VALUE' : 'VALUE'
            }
    
    # Sets the capacity of the line to the average of max and max counter flows.
    df['VALUE'] = (df['max_flow'] + df['max_counter_flow']) / 2
    
    df = df[list(cols.keys())]
    
    df = df.rename(columns = cols)
    
    return df

def format_gtd_planned(df):
    cols = {'pathway' : 'TECHNOLOGY',
            'from_region' : 'From',
            'to_region' : 'To', 
            'year_planned' : 'YEAR',
            'VALUE' : 'VALUE'
            }
    
    # Sets the capacity of the line to the average of max and max counter flows.
    df['VALUE'] = (df['max_flow'] + df['max_counter_flow']) / 2
    
    df = df[list(cols.keys())]
    
    df = df.rename(columns = cols)
    
    return df

def correct_gtd_data(df_exist, df_plan, region_mapping_dict, 
                     CUSTOM_TRN_BA_DICT_FROM, CUSTOM_TRN_BA_DICT_TO,
                     CUSTOM_TRN_BA_MISSING):
    
    # Create df with unique transmission technologies from GTD.
    tech_exist = get_unique_technologies(df_exist)
    tech_plan = get_unique_technologies(df_plan)
    tech_trn_list = list(set(tech_exist + tech_plan)) 
    df = pd.DataFrame(tech_trn_list, columns = ['TECHNOLOGY_gtd'])
 
    # Map OG regions to GTD technologies.
    df['From'], df['To'] = df['TECHNOLOGY_gtd'].str[3:8], df['TECHNOLOGY_gtd'].str[8:13]
    df['From'], df['To'] = df['From'].map(region_mapping_dict
                                          ), df['To'].map(region_mapping_dict)
    
    """For USA GTD regions that are larger than certain OG regions (MISO, PJM, SWPP),
    as well as for the Nei Mongol region in China, assign custom OG regions as a 
    best estimate for the most important transmission pathway. This pathway will have
    all residual capacity assigned to it."""
    df.loc[df['From'] == 'DUPLICATE', 
           'From'] = df['TECHNOLOGY_gtd'].map(CUSTOM_TRN_BA_DICT_FROM)
    
    df.loc[df['To'] == 'DUPLICATE', 
           'To'] = df['TECHNOLOGY_gtd'].map(CUSTOM_TRN_BA_DICT_TO)
    
    # Sort 'From' and 'To' columns alphabetically for custom entries. 
    for idx,row in df.loc[df['To'] < df['From']].iterrows():
        df.loc[df['TECHNOLOGY_gtd'] == row['TECHNOLOGY_gtd'], 
               ['From', 'To']] = row['To'], row['From']

    df.set_index('TECHNOLOGY_gtd', inplace = True)
    
    # Create corrected TECHNOLOGY entry
    df['TECHNOLOGY'] = 'TRN' + df['From'] + df['To']

    # Update original GTD df's with corrected entries
    df_exist.set_index('TECHNOLOGY', inplace = True)
    df_exist.update(df)
    df_exist = df_exist.merge(df[['TECHNOLOGY']], left_index = True, 
                              right_index = True).reset_index(drop = True)
    df_exist = df_exist[df_exist['From'] != df_exist['To']]
    
    df_plan.set_index('TECHNOLOGY', inplace = True)
    df_plan.update(df)
    df_plan = df_plan.merge(df[['TECHNOLOGY']], left_index = True, 
                            right_index = True).reset_index(drop = True)   
    df_plan = df_plan[df_plan['From'] != df_plan['To']]
    
    # Group corrected GTD entries
    df_exist = df_exist.groupby(['TECHNOLOGY', 'From', 'To']
                                ).sum('VALUE').reset_index()
    
    df_plan = df_plan.groupby(['TECHNOLOGY', 'From', 'To', 'YEAR']
                                ).sum('VALUE').reset_index()
    
    """Add missing transmission pathways following the approach to disaggregate large GTD
    regions (e.g. MISO, PJM, SWPP). Transmission pathways between OG regions that fall 
    within a single GTD region (e.g. USASN and USASS in MISO) are not automatically 
    pulled from the GTD datasets so need to be manually added. Furthermore, pathways
    within the GTD datasets that are not selected as the pathways to which residual
    capacity should be allocated (defined in the CUSTOM_TRN_BA_DICT_FROM and
    CUSTOM_TRN_BA_DICT_TO constants) also need to be manually added."""
    df_missing = pd.DataFrame(CUSTOM_TRN_BA_MISSING, columns = ['TECHNOLOGY'])
    df_missing['From'] = df_missing['TECHNOLOGY'].str[3:8]
    df_missing['To'] = df_missing['TECHNOLOGY'].str[8:13]
    df_missing['YEAR'] = '-'
    df_missing['VALUE'] = 0

    df_exist = pd.concat([df_exist, df_missing], join = 'inner').reset_index(drop = True)
    df_plan = pd.concat([df_plan, df_missing]).reset_index(drop = True)
    
    return df_exist, df_plan

def calculate_transmission_distances(df_exist_corrected, df_plan_corrected, 
                                     centerpoints_dict):
    
    # Create df with unique transmission technologies from corrected GTD data.
    df = pd.concat([df_exist_corrected[['TECHNOLOGY', 'From', 'To']], 
                    df_plan_corrected[['TECHNOLOGY', 'From', 'To']]]).drop_duplicates()
    
    # Add coordinates to region entries.
    for row in centerpoints_dict:
      df.loc[df['From'] == row['region'], 'From_lat'] = row['lat']
      df.loc[df['From'] == row['region'], 'From_long'] = row['long']
      
      df.loc[df['To'] == row['region'], 'To_lat'] = row['lat']
      df.loc[df['To'] == row['region'], 'To_long'] = row['long']
    
    # Calculate respetive transmission distances.
    df['distance'] = df.apply(lambda x: geopy.distance.geodesic((x['From_lat'], x['From_long']), 
                                                          (x['To_lat'], x['To_long'])), axis=1)
    
    df['distance'] = df['distance'].astype(str).str.replace(' km', '').astype(float).round(0)
    
    return df

def set_break_even_distance(trn_param):
    
    # Set base non-distance dependent CAPEX costs (converter pair).
    hvac_cost = trn_param['HVAC'][1]
    hvdc_cost = trn_param['HVDC'][1]
    n = 1
    
    """Iterate over distance dependent CAPEX costs (line) until the
    break-even distance is found after which HVDC is cheaper"""
    while hvdc_cost > hvac_cost:
        hvac_cost = hvac_cost + trn_param['HVAC'][0]
        hvdc_cost = hvdc_cost + trn_param['HVDC'][0]
    
        n = n + 1
        
    return n

def set_transmission_tech_groups(df_exist_corrected, df_plan_corrected, 
                                 centerpoints_dict, trn_param, SUBSEA_LINES):
    
    # Set the break-even distance.
    break_even_dist = set_break_even_distance(trn_param)
    
    # Set the transmission distances.
    df = calculate_transmission_distances(df_exist_corrected, df_plan_corrected, 
                                          centerpoints_dict)
    
    """Compare the break-even distance with transmission distances by pathway and set
    the respective technology group per individual transmission pathway."""
    df.loc[df['distance'] >= break_even_dist, 
           'tech_group'] = 'HVDC'
    
    df.loc[df['distance'] < break_even_dist, 
           'tech_group'] = 'HVAC'
    
    df.loc[df['TECHNOLOGY'].isin(SUBSEA_LINES) , 
           'tech_group'] = 'HVDC_subsea'
    
    tech_group_dict = dict(zip(df['TECHNOLOGY'], 
                               df['tech_group']))
    
    return tech_group_dict