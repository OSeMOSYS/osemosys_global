# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 14:50:11 2022

@author: mbrinkerink
"""
import urllib
import os
import yaml
import pandas as pd

yaml_file = open("config.yaml")
parsed_yaml_file = yaml.load(yaml_file, Loader=yaml.FullLoader)

input_dir = parsed_yaml_file.get('inputDir')

# %% COPY PASTE BELOW AND INSERT IN OPG_powerplant_data.py

## Checks whether PLEXOS-World/MESSAGEix-GLOBIOM soft-link model data needs to be retrieved from the PLEXOS-World Harvard Dataverse.
try:
    Open = pd.read_excel(os.path.join(input_dir,
                             "PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx"))

except IOError:
    
    url = 'https://dvn-cloud.s3.amazonaws.com/10.7910/DVN/O6ICJP/17f2cb5b2ff-174d6e018d76?response-content-disposition=attachment%3B%20filename%2A%3DUTF-8%27%27PLEXOS-World%2520model%2520MESSAGEix%2520-%2520GLOBIOM%2520Soft-Link.xlsx&response-content-type=application%2Fvnd.openxmlformats-officedocument.spreadsheetml.sheet&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20220225T173557Z&X-Amz-SignedHeaders=host&X-Amz-Expires=3600&X-Amz-Credential=AKIAIEJ3NV7UYCSRJC7A%2F20220225%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=b5aee3a0082723fbce0db78512e5cfdae9fc267ac568a46d103842fdf79d1cf2'

    urllib.request.urlretrieve(url, os.path.join(input_dir, 'PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx'))
    
df_reslimit = pd.read_excel(os.path.join(input_dir, "PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx"), sheet_name = "Properties")

dict_reslimit = {'Hydro' : 'HYD',
                 'Solar|CSP' : 'CSP',
                 'Solar|PV' : 'SPV',
                 'Wind|Onshore' : 'WON',
                 'Wind|Offshore' : 'WOF'}

df_reslimit_units = df_reslimit.loc[(df_reslimit['child_object'].str.contains('|'.join(dict_reslimit.keys()))) & 
                              (df_reslimit['property'] == 'Max Units Built') &
                              (df_reslimit['scenario'].str.contains('Base')) &
                              (df_reslimit['child_class'] == 'Generator')].set_index('child_object')

df_reslimit_capacity = df_reslimit.loc[(df_reslimit['child_object'].str.contains('|'.join(dict_reslimit.keys()))) & 
                              (df_reslimit['property'] == 'Max Capacity') &
                              (df_reslimit['child_class'] == 'Generator')].set_index('child_object')

df_reslimit_final = pd.DataFrame(df_reslimit_capacity['value'] * df_reslimit_units['value'] / 1000).rename(columns = {'value' : 'VALUE'})
df_reslimit_final['node'], df_reslimit_final['powerplant']  = df_reslimit_final.index.str.rsplit('|',1).str[1], df_reslimit_final.index.str.rsplit('|',1).str[0]
df_reslimit_final['powerplant'] = df_reslimit_final['powerplant'].map(dict_reslimit)

df_reslimit_final.loc[df_reslimit_final['node'].str.len() <= 6, 
             'node_code'] = (df_reslimit_final['node'].
                             str.split('-').
                             str[1:].
                             str.join("") +
                             'XX')
df_reslimit_final.loc[df_reslimit_final['node'].str.len() > 6, 
             'node_code'] = (df_reslimit_final['node'].
                             str.split('-').
                             str[1:].
                             str.join("")
                            )

df_reslimit_final['TECHNOLOGY'] = 'PWR' + df_reslimit_final['powerplant'] + df_reslimit_final['node_code'] + '01'
df_reslimit_final = df_reslimit_final[['TECHNOLOGY', 'VALUE']].reset_index(drop = True)