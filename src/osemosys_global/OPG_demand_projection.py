#!/usr/bin/env python
# coding: utf-8

# # Script for electricity demand projection and downscaling

# In[1]:


import pandas as pd
import world_bank_data as wb
import numpy as np
import matplotlib.pyplot as plt
import urllib
import xlsxwriter
import os
from sklearn.linear_model import LinearRegression
import yaml


# ## Input data and projection boundaries

# ### Imports relevant input data

# In[2]:

#Read in information from YAML file
yaml_file = open("config.yaml")
parsed_yaml_file = yaml.load(yaml_file, Loader=yaml.FullLoader)

input_dir = parsed_yaml_file.get('inputDir')
output_dir = parsed_yaml_file.get('outputDir') + 'data/'                

#Checks whether PLEXOS-World 2015 data needs to be retrieved from the PLEXOS-World Harvard Dataverse.

try:
    Open = open(os.path.join(input_dir,'PLEXOS_World_2015_Gold_V1.1.xlsx'))
    
except IOError:
    urllib.request.urlretrieve("https://dataverse.harvard.edu/api/access/datafile/4008393?format=original&gbrecs=true" , 
                               os.path.join(input_dir,"PLEXOS_World_2015_Gold_V1.1.xlsx"))
    
    Open = open(os.path.join(input_dir,"PLEXOS_World_2015_Gold_V1.1.xlsx"))

finally:
    Open.close()

Import_memberships = pd.read_excel(os.path.join(input_dir,"PLEXOS_World_2015_Gold_V1.1.xlsx") , sheet_name = "Memberships")

#Imports SSP GDPppp and Population projections (https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=30)
Import_iamc_db_GDPppp_Countries = pd.read_excel(os.path.join(input_dir,'iamc_db_GDPppp_Countries.xlsx'))
Import_iamc_db_POP_Countries = pd.read_excel(os.path.join(input_dir,'iamc_db_POP_Countries.xlsx'))
Import_iamc_db_URB_Countries = pd.read_excel(os.path.join(input_dir,'iamc_db_URB_Countries.xlsx'))

#Imports custom GDPppp and Population projections for countries not included in the SSP datasets.
Import_POP_Missing = pd.read_excel(os.path.join(input_dir,'iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx') , 
                                   sheet_name = 'POP').set_index('Region')

Import_GDP_Missing = pd.read_excel(os.path.join(input_dir,'iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx') , 
                                   sheet_name = 'GDP|PPP').set_index('Region')

Import_URB_Missing = pd.read_excel(os.path.join(input_dir,'iamc_db_POP_GDPppp_URB_Countries_Missing.xlsx') , 
                                   sheet_name = 'URB').set_index('Region')

#Imports T&D losses projections (https://www.sciencedirect.com/science/article/pii/S0142061518335075?via%3Dihub)
Import_Incl_Losses = pd.read_excel(os.path.join(input_dir,'T&D Losses.xlsx'))

Import_iamc_db_GDPppp_Countries.head(2)


# ### Set boundaries for the regression

# In[3]:


#Sets the spatial resolution for the regession, right now can only be applied to continental level by setting 'child_object'. 
#Note that not all countries have historical data available so country-level regression can't be applied for all.
Spatial_Resolution = 'child_object'

#Include urbanization as part of the regression? 'Yes' or 'No'
Urbanization = 'Yes'

#Set which SSP data and sources are to be used for the country-level GDP|PPP and Population projects.
Pathway = 'SSP2'
POP_Countries_Source = 'IIASA-WiC POP' #Options are 'IIASA-WiC POP' and 'NCAR'
GDPppp_Countries_Source = 'OECD Env-Growth' #Options are 'IIASA GDP' and 'OECD Env-Growth'. 'OECD Env-Growth' has more entries.
URB_Countries_Source = 'NCAR' #'NCAR' is the only option.

#Set modelling years
Years_List_5 = np.arange(2010, 2105, 5).tolist()
Years_List = np.arange(2010, 2101).tolist()
Years_List_Interpolation = list(set(Years_List) - set(Years_List_5))


# In[4]:


Spatial_Mapping = Import_memberships.copy()
Spatial_Mapping = Spatial_Mapping[Spatial_Mapping['collection'].str.contains('Region')]
Spatial_Mapping['Country'] = Spatial_Mapping.parent_object.str.split('-' , expand = True)[1]

Spatial_Mapping_Country = Spatial_Mapping.copy().set_index('Country')
Spatial_Mapping_Country = Spatial_Mapping_Country.loc[~Spatial_Mapping_Country.index.duplicated(keep='first')]

Spatial_Mapping_Node = Spatial_Mapping.copy().set_index('parent_object')
Spatial_Mapping_Node = Spatial_Mapping_Node.loc[~Spatial_Mapping_Node.index.duplicated(keep='first')]
               
Spatial_Mapping_Node.head(1)


# ### Retrieves PLEXOS-World 2015 hourly demand data incl. T&D losses for all nodes as baseline value for the demand forecasting
# Used to be able to disaggregate regional electricity demand to the nodal level as well as calculate relative peak demand per node.

# In[6]:


#Checks whether PLEXOS-World 2015 data needs to be retrieved from the PLEXOS-World Harvard Dataverse.
try:
    Open = open(os.path.join(input_dir,'All_Demand_UTC_2015.csv'))
    
    Import_Hourly_Demand_2015 = pd.read_csv(os.path.join(input_dir,'All_Demand_UTC_2015.csv') , encoding='latin-1')
    
except IOError:
    urllib.request.urlretrieve ('https://dataverse.harvard.edu/api/access/datafile/3985039?format=original&gbrecs=true', 
                                os.path.join(input_dir,'All_Demand_UTC_2015.csv'))
    
    Import_Hourly_Demand_2015 = pd.read_csv(os.path.join(input_dir,'All_Demand_UTC_2015.csv') , encoding='latin-1')
    
Import_Hourly_Demand_2015.head(2)


# ### Determines relative 2015 share of demand per sub-country node

# In[10]:


#Sums the hourly demand as retrieved from the PLEXOS-World dataset to year total (in MWh) and drops all hourly values.
Demand_2015_Raw = Import_Hourly_Demand_2015.drop(columns = ['Datetime'])

Demand_2015_Raw.loc['Node_Demand_2015'] = Demand_2015_Raw.sum()

Demand_2015_Raw = Demand_2015_Raw.iloc[8760:]

#Transposes the dataframe and uses the original headers as column entry. 
Demand_2015_Raw = Demand_2015_Raw.transpose().reset_index().rename(columns = {'index' : 
                                                                              'PLEXOS_Nodes'})

#Adds country entry to dataframe (e.g. NA-USA-CA Node gets column entry NA-USA country)
Demand_2015_Raw.insert(loc = 0 , 
                       column = 'Country' , 
                       value  = Demand_2015_Raw.PLEXOS_Nodes.str.split('-' , expand = True)[1])

Demand_2015_Raw.insert(loc = 1 , 
                       column = 
                       'PLEXOS_Countries' , 
                       value = Demand_2015_Raw['PLEXOS_Nodes'].str[:6])

#Creates a dataframe excluding all sub-country nodes
Country_Demand_2015 = Demand_2015_Raw[Demand_2015_Raw['PLEXOS_Countries'] == 
                                      Demand_2015_Raw['PLEXOS_Nodes']
                                     ].drop(columns = ['PLEXOS_Nodes']
                                           ).rename(columns = {'Node_Demand_2015' : 
                                                            'Country_Demand_2015'})

#Adds country-level 2015 demand in column adjacent to sub-country level 2015 demand and calculates relative share per country.
Node_Demand_2015 = pd.merge(Demand_2015_Raw , 
                            Country_Demand_2015[['PLEXOS_Countries' , 
                                                 'Country_Demand_2015']] , 
                            on = 'PLEXOS_Countries' , 
                            how ='left')

Node_Demand_2015['Share_%_Country_Demand'] = Node_Demand_2015['Node_Demand_2015'
                                                             ] / Node_Demand_2015['Country_Demand_2015']

Node_Demand_2015.iloc[49:50]


# ## Historical relationships

# ### Creates historic relationships based on World Bank Data

# In[11]:


#Extracts historical GDPppp per capita (constant 2017 international $) from the World Bank API
Country_GDPppp_WB = wb.get_series('NY.GDP.PCAP.PP.KD', 
                                  date = '1980:2014' , 
                                  id_or_value = 'id')

Country_GDPppp_WB = Country_GDPppp_WB.reset_index().rename(columns = 
                                                           {'NY.GDP.PCAP.PP.KD' : 
                                                            'WB_GDPppp'}
                                                          ).set_index('Country')

#Extracts Electricity consumption per capita in kWh from the World Bank API. Data available for up till 2014
Country_Elec_WB = wb.get_series('EG.USE.ELEC.KH.PC', 
                                date = '1980:2014' , 
                                id_or_value = 'id')

Country_Elec_WB = Country_Elec_WB.reset_index().rename(columns = 
                                                       {'EG.USE.ELEC.KH.PC' : 
                                                        'WB_Elec'}
                                                      ).set_index('Country')

#Extracts Urban population (% of total population) from the World Bank API.
Country_Urb_WB = wb.get_series('SP.URB.TOTL.IN.ZS', 
                               date = '1980:2014' , 
                               id_or_value = 'id')

Country_Urb_WB = Country_Urb_WB.reset_index().rename(columns = 
                                                     {'SP.URB.TOTL.IN.ZS' : 
                                                      'WB_Urb'}
                                                    ).set_index('Country')

Country_Urb_WB.tail(1)


# ### Applies linear regression

# In[12]:


#Merges the relevant dataframes
Country_Regression_WB = pd.merge(Country_GDPppp_WB[['Year' , 
                                                     'WB_GDPppp']] , 
                                  Country_Elec_WB[['Year' , 
                                                   'WB_Elec']] , 
                                  left_on = ['Year' , 
                                             'Country'] , 
                                  right_on = ['Year' , 
                                              'Country'])

if Urbanization == 'Yes':
    Country_Regression_WB = pd.merge(Country_Regression_WB[['Year' , 
                                                            'WB_GDPppp' , 
                                                            'WB_Elec']] , 
                                     Country_Urb_WB[['Year' , 
                                                     'WB_Urb']] , 
                                  left_on = ['Year' , 
                                             'Country'] , 
                                  right_on = ['Year' , 
                                              'Country'])

#Drops all entries that don't have an inner match
Country_Regression_WB = Country_Regression_WB.dropna()

Country_Regression_WB = pd.merge(Country_Regression_WB, 
                                  Spatial_Mapping_Country[['parent_object' , 
                                                           'child_object']] , 
                                  left_index = True , 
                                  right_index = True)

Country_Regression_WB = Country_Regression_WB.set_index(Spatial_Resolution)

#Groups the entries by <Spatial_Resolution> and calculates the regional linear fit based on all historical values 
sklearn_lr = LinearRegression()
Country_Regression_WB_Grouped = pd.DataFrame()

for x in Country_Regression_WB.index.unique():
    
    Country_Regression_WB_Temp = Country_Regression_WB.copy().loc[x]
    
    #If Urbanization is included linear regression occurs with multiple independent variables (GDPppp and % Urban population)
    #for the dependent variable (Electricity demand). 
    if Urbanization == 'Yes':
        
        sklearn_lr.fit(Country_Regression_WB_Temp[['WB_GDPppp' , 
                                                   'WB_Urb']] , 
                       Country_Regression_WB_Temp['WB_Elec'])
        
        Country_Regression_WB_Temp['intercept'
                                  ] ,  Country_Regression_WB_Temp['coef_GDPppp'
                                                                 ] ,  Country_Regression_WB_Temp['coef_Urb'
                                                                                                ] = [sklearn_lr.intercept_ , 
                                                                                                float(sklearn_lr.coef_[0]) , 
                                                                                               float(sklearn_lr.coef_[1])]
               
        Country_Regression_WB_Temp['R2_GDPppp_Urb/Elec'] = sklearn_lr.score(Country_Regression_WB_Temp[['WB_GDPppp' , 
                                                   'WB_Urb']] , Country_Regression_WB_Temp['WB_Elec'])
        
        Country_Regression_WB_Grouped = Country_Regression_WB_Grouped.append(Country_Regression_WB_Temp)
        
    #If Urbanization is not included linear regression occurs with single independent variables (GDPppp) for the dependent
    #variable (Electricity demand).    
    else:
        
        sklearn_lr.fit(Country_Regression_WB_Temp[['WB_GDPppp']], 
                       Country_Regression_WB_Temp['WB_Elec'])
              
        Country_Regression_WB_Temp['intercept'
                                  ] ,  Country_Regression_WB_Temp['coef_GDPppp'] = [sklearn_lr.intercept_ , 
                                                                                    float(sklearn_lr.coef_)]
        
        Country_Regression_WB_Temp['R2_GDPppp/Elec'] = sklearn_lr.score(Country_Regression_WB_Temp[['WB_GDPppp']] , 
                                                                   Country_Regression_WB_Temp['WB_Elec'])
        
        Country_Regression_WB_Grouped = Country_Regression_WB_Grouped.append(Country_Regression_WB_Temp)

Country_Regression_WB_Grouped.head(2)


# In[13]:


for a in Country_Regression_WB_Grouped.index.unique():
    Country_Regression_WB_plot = Country_Regression_WB_Grouped.loc[a]
    b = Country_Regression_WB_plot['WB_GDPppp']
    c = Country_Regression_WB_plot['WB_Elec']
    x = Country_Regression_WB_plot.loc[a , 'WB_GDPppp']
    m = Country_Regression_WB_plot.loc[a , 'coef_GDPppp']
    z = Country_Regression_WB_plot.loc[a , 'intercept']
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(b , c , color = 'blue' , alpha = 0.5 , label = 'GDPppp per capita')
    
    if Urbanization == 'Yes':
        
        w2 = Country_Regression_WB_plot.loc[a , 'R2_GDPppp_Urb/Elec'].unique()
        d = Country_Regression_WB_plot['WB_Urb']
        ax2 = ax1.twiny()
        ax2.scatter(d , c , color = 'red' , alpha = 0.5 , label = 'Urban population') 
        ax2.set_xlabel('Urban population (% of total population)')
        ax2.legend(loc = 'upper right')
        plt.title(f'{a} R2 = {w2}')
        
    else:
        ax1.plot(x , m * x + z , color = 'black' , alpha = 0.5)
        w = Country_Regression_WB_plot.loc[a , 'R2_GDPppp/Elec'].unique()
        plt.title(f'{a} R2 = {w}')
    
    ax1.legend(loc='upper left')
    ax1.set_xlabel('GDPppp per capita (constant 2017 international $)')
    ax1.set_ylabel('Electricity demand per capita (kWh)')
    plt.rcParams["figure.figsize"] = [8,4]
    
    # Rather than plt.show() we need to make this save...
    # plt.show()
    plt.savefig(os.path.join(output_dir,r'../figs',f'{a} R2 = {w2}'+'.jpg'))


# ## Country-level Projections

# ### Creates dataframe with SSP specific population projections

# In[14]:


#Filters relevant data based on earlier given entries (SSP Pathway and source for POP)
Country_POP_SSP = Import_iamc_db_POP_Countries.loc[(Import_iamc_db_POP_Countries['Model'] == 
                                                    POP_Countries_Source) 
                                                   & (Import_iamc_db_POP_Countries['Scenario'] == 
                                                      Pathway)
                                                  ].set_index('Region')

#Checks whether pathway data is available for all included countries. In case pathway data is not available for all countries 
#it checks wheter custom data is provided in 'iamc_db_POP_POP_Countries_Missing.xlsx'. If not, an error POPs up indicating 
#for which countries data is missing. Data has to be manually added to project demand for all countries.
Country_POP_SSP_Missing = Spatial_Mapping_Country[(~Spatial_Mapping_Country.index.isin
                                                   (Country_POP_SSP.index))]

for x in Country_POP_SSP_Missing.index:
    a = x in Import_POP_Missing.index
    if a == False:
        raise SystemExit(f'Country data for {x} is not available in the SSP or custom dataset! Add country-level data.')
    else: 
        print(f'Country data for {x} included from custom dataset.')

#Filters data for relevant SSP
Import_POP_Missing = Import_POP_Missing.loc[(Import_POP_Missing['Scenario'] == Pathway)]

#Appends both dataframes
Country_POP_SSP = Country_POP_SSP.append(Import_POP_Missing, ignore_index = False)

#Filters data for relevant to be modelled countries
Country_POP_SSP = pd.merge(Spatial_Mapping_Country[['child_object']] , 
                           Country_POP_SSP , 
                           left_index = True , 
                           right_index = True , 
                           how = 'inner')

Country_POP_SSP.tail(1)


# ### Creates dataframe with SSP specific GDP|PPP projections

# In[15]:


#Filters relevant data based on earlier given entries (SSP Pathway and source for GDPppp)
Country_GDPppp_SSP = Import_iamc_db_GDPppp_Countries.loc[(Import_iamc_db_GDPppp_Countries['Model'] == 
                                                          GDPppp_Countries_Source) & 
                                                         (Import_iamc_db_GDPppp_Countries['Scenario'] == 
                                                          Pathway)
                                                        ].set_index('Region')

#Checks whether pathway data is available for all included countries. In case pathway data is not available for all countries 
#it checks wheter custom data is provided in 'iamc_db_GDPppp_GDPppp_Countries_Missing.xlsx'. If not, an error GDPppps up 
#indicating for which countries data is missing. Data has to be manually added to project demand for all countries.
Country_GDPppp_SSP_Missing = Spatial_Mapping_Country[(~Spatial_Mapping_Country.index.isin(Country_GDPppp_SSP.index))]

for x in Country_GDPppp_SSP_Missing.index:
    a = x in Import_GDP_Missing.index
    if a == False:
        raise SystemExit(f'Country data for {x} is not available in the SSP or custom dataset! Add country-level data.')
    else: 
        print(f'Country data for {x} included from custom dataset.')

#Filters data for relevant SSP
Import_GDP_Missing = Import_GDP_Missing.loc[(Import_GDP_Missing['Scenario'] == Pathway)]

#Appends both dataframes
Country_GDPppp_SSP = Country_GDPppp_SSP.append(Import_GDP_Missing, ignore_index = False)

#Filters data for relevant to be modelled countries
Country_GDPppp_SSP = pd.merge(Spatial_Mapping_Country[['child_object']] , 
                              Country_GDPppp_SSP , 
                              left_index = True , 
                              right_index = True , 
                              how = 'inner')

Country_GDPppp_SSP.head(1)


# ### Creates dataframe with SSP specific urban population projections

# In[16]:


#Filters relevant data based on earlier given entries (SSP Pathway and source for URB)
Country_URB_SSP = Import_iamc_db_URB_Countries.loc[(Import_iamc_db_URB_Countries['Model'] == 
                                                          URB_Countries_Source) & 
                                                         (Import_iamc_db_URB_Countries['Scenario'] == 
                                                          Pathway)
                                                        ].set_index('Region')

#Checks whether pathway data is available for all included countries. In case pathway data is not available for all countries 
#it checks wheter custom data is provided in 'iamc_db_URB_URB_Countries_Missing.xlsx'. If not, an error URBs up 
#indicating for which countries data is missing. Data has to be manually added to project demand for all countries.
Country_URB_SSP_Missing = Spatial_Mapping_Country[(~Spatial_Mapping_Country.index.isin(Country_URB_SSP.index))]

for x in Country_URB_SSP_Missing.index:
    a = x in Import_URB_Missing.index
    if a == False:
        raise SystemExit(f'Country data for {x} is not available in the SSP or custom dataset! Add country-level data.')
    else: 
        print(f'Country data for {x} included from custom dataset.')

#Filters data for relevant SSP
Import_URB_Missing = Import_URB_Missing.loc[(Import_URB_Missing['Scenario'] == Pathway)]

#Appends both dataframes
Country_URB_SSP = Country_URB_SSP.append(Import_URB_Missing, ignore_index = False)

#Filters data for relevant to be modelled countries
Country_URB_SSP = pd.merge(Spatial_Mapping_Country[['child_object']] , 
                              Country_URB_SSP , 
                              left_index = True , 
                              right_index = True , 
                              how = 'inner')

Country_URB_SSP.head(1)


# ### Projects electricity demand by making use of historic relationships and SSP specific Population, GDP|PPP and optionally urbanization projections

# In[17]:


#Creates base dataframe for demand projections with the required coefficients
Base_df_projections = Country_GDPppp_SSP[['child_object' , 
                                          'Scenario']]

Regression_Coefficients = Country_Regression_WB_Grouped.loc[~Country_Regression_WB_Grouped.index.duplicated(keep = 'first')]

if Urbanization == 'Yes':
    w = ['coef_GDPppp' , 
         'coef_Urb' , 
         'intercept']
else:
    w = ['coef_GDPppp' , 
         'intercept']
    
#Divides the country level GDP|PPP data (converted from billions to millions) with the population (millions) to get GDP|PPP pp.
Country_GDPppp_pp_SSP = pd.merge(Base_df_projections , 
                           Regression_Coefficients[w] , 
                           left_on = 'child_object' , 
                           right_index = True , 
                           how = 'left')

for x in Years_List_5:
    Country_GDPppp_pp_SSP[x] = (Country_GDPppp_SSP[x] * 1000) / Country_POP_SSP[x]
    
if Urbanization == 'Yes':
    y = ['child_object' , 
         'Scenario' , 
         'coef_GDPppp' ,
         'coef_Urb' ,
         'intercept']
else:
    y = ['child_object' , 
             'Scenario' , 
             'coef_GDPppp' ,
             'intercept']
    
#Applies the country-level GDP|PPP pp values to a linear regression to project country-level electricity demand pp.
Country_Demand_projected_SSP = Country_GDPppp_pp_SSP.copy([[y]])

Country_Demand_projected_SSP['Variable'] = 'Demand|projected|pp'

if Urbanization == 'Yes':
    for z in Years_List_5:
         Country_Demand_projected_SSP[z] = (Country_Demand_projected_SSP['coef_GDPppp'] * 
                                            Country_GDPppp_pp_SSP[z] + 
                                            Country_Demand_projected_SSP['coef_Urb'] * 
                                            Country_URB_SSP[z] + 
                                            Country_Demand_projected_SSP['intercept'])
else:
    for z in Years_List_5:
         Country_Demand_projected_SSP[z] = (Country_Demand_projected_SSP['coef_GDPppp'] * 
                                            Country_GDPppp_pp_SSP[z] + 
                                            Country_Demand_projected_SSP['intercept'])
        
Country_Demand_projected_SSP.head(1)


# In[18]:


for a in Spatial_Mapping_Country['child_object'].unique():
    Country_Regression_WB_Grouped_plot = Country_Regression_WB_Grouped.loc[a]
    
    Country_GDPppp_pp_SSP_plot = Country_GDPppp_pp_SSP.loc[(Country_GDPppp_pp_SSP['child_object'] == a)]
    
    Country_Demand_projected_SSP_plot = Country_Demand_projected_SSP.loc[(Country_Demand_projected_SSP['child_object'
                                                                                                         ] == a)]
    x = Country_Regression_WB_Grouped_plot['WB_GDPppp']
    x2 = Country_GDPppp_pp_SSP_plot[2035]
    x3 = Country_GDPppp_pp_SSP_plot[2050]
    x4 = Country_GDPppp_pp_SSP_plot[2100]
    y = Country_Regression_WB_Grouped_plot['WB_Elec']
    y2 = Country_Demand_projected_SSP_plot[2035]
    y3 = Country_Demand_projected_SSP_plot[2050]
    y4 = Country_Demand_projected_SSP_plot[2100]
    m = Country_Regression_WB_Grouped_plot['coef_GDPppp'].unique()
    b = Country_Regression_WB_Grouped_plot['intercept'].unique()
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(x , y , color = 'grey' , alpha = 0.5 , label = '1980-2014 (World Bank)')
    ax1.scatter(x2 , y2 , color = 'green' , alpha = 0.5 , label = '2035')
    ax1.scatter(x3 , y3 , color = 'blue' , alpha = 0.5 , label = '2050')
    ax1.scatter(x4 , y4 , color = 'red' , alpha = 0.5 , label = '2100')
    ax1.legend(loc='upper left')
    
    if Urbanization == 'No':
        plt.plot(x , m * x + b , color = 'black' , alpha = 0.5)
        plt.plot(x2 , m * x2 + b , color = 'black' , alpha = 0.5)
        plt.plot(x3 , m * x3 + b , color = 'black' , alpha = 0.5)
        plt.plot(x4 , m * x4 + b , color = 'black' , alpha = 0.5)
        
    plt.title(f'Demand projection {a}')
    ax1.set_xlabel('GDPppp per capita (constant 2017 international $)')
    ax1.set_ylabel('Electricity demand per capita (kWh)')
    plt.rcParams["figure.figsize"] = [8,4]
    
    # Rather than plt.show() we need to make this save...
    # plt.show()
    plt.savefig(os.path.join(output_dir,r'../figs',f'Demand projection {a}'+'.jpg'))


# ### Aggregates projected demand per person to full country-level

# In[19]:


#Multiplies the country-level projected demand pp (in kWh) with the total population (in millions) to get country-level 
#total projected demand (in GWh).
Country_Demand_projected_SSP_Aggregated = pd.DataFrame()
for x in Years_List_5:
     Country_Demand_projected_SSP_Aggregated[x] = (Country_Demand_projected_SSP[x] * 
                                                   Country_POP_SSP[x]).round(2)
        
Country_Demand_projected_SSP_Aggregated.head(1)


# ### Adds transmission and distribution losses to country-level demand
# Explicit modelling of domestic transmission and distribution is not incorporated in PLEXOS-World. Country-level T&D losses per 5-year interval are added to the projected electricity demand based on Sadovskaia et al., 2019; https://doi.org/10.1016/j.ijepes.2018.11.012. Study includes data for up till 2050. 5-year intervals after that are manually added (Maarten Brinkerink) and values kept equal compared to 2050.

# In[20]:


#Checks whether T&D is available for all included countries. An assertion error pops up in case data is missing. 
#In this case either data needs to be added to 'T&D Losses.xlsx'.
Losses_Missing = Spatial_Mapping_Country[(~Spatial_Mapping_Country.index.isin(Import_Incl_Losses.Country))]

assert Losses_Missing.empty , print(Losses_Missing)

Country_Demand_Incl_Losses = Import_Incl_Losses.copy().set_index('Country')

Country_Demand_Incl_Losses = Country_Demand_Incl_Losses[Country_Demand_Incl_Losses.index.isin
                                                        (Spatial_Mapping_Country.index)
                                                       ].reindex(Spatial_Mapping_Country.index)

Country_Demand_projected_SSP_Incl_Losses_Raw = pd.DataFrame()

#Add T&D losses to the projected country-level demand.
for z in Years_List_5:
     Country_Demand_projected_SSP_Incl_Losses_Raw[z] = (Country_Demand_projected_SSP_Aggregated[z] * 
                                                    Country_Demand_Incl_Losses[z] / 100 + 
                                                    Country_Demand_projected_SSP_Aggregated[z]
                                                   ).round(2).astype(float)

Country_Demand_projected_SSP_Incl_Losses_Raw.tail(1)


# ### Constraints the forecasted final demand to 2015 baseline values as minimum
# In case of linear regression, smaller countries with signficantly lower projected independent variables (GDP, Urbanization) compared to the regional average can lead to very low and often negative projected demand values (e.g. EU-KOS). Hence, a comparison is being made to the 2015 baseline demand values with the assumption that a decline in electricity demand is not realistic (note: as of now no decoupling of GDP growth and energy demand reduction has been assumed).

# In[21]:


Country_Demand_projected_SSP_Incl_Losses = pd.DataFrame()

for x in Country_Demand_projected_SSP_Incl_Losses_Raw.index.unique():
    Country_Demand_projected_SSP_Incl_Losses_x = Country_Demand_projected_SSP_Incl_Losses_Raw.loc[x]
    
    Node_Demand_2015_x = Node_Demand_2015[Node_Demand_2015['Country'].str.contains(x)]

    Country_Demand_projected_SSP_Incl_Losses_x = Country_Demand_projected_SSP_Incl_Losses_x.clip(
        (Node_Demand_2015_x.iloc[0]['Country_Demand_2015'] / 1000) , )
    
    Country_Demand_projected_SSP_Incl_Losses = Country_Demand_projected_SSP_Incl_Losses.append(
        Country_Demand_projected_SSP_Incl_Losses_x)

Country_Demand_projected_SSP_Incl_Losses.tail(1)


# ## Downscaling from country-level to nodal-level

# ### Uses relative 2015 share in demand per sub-country node to downscale country-level scenario specific demand

# In[22]:


#Drops the non-required columns and merges the 2015 and projected demand dataframes
Node_Demand_SSP_projected_Incl_Losses  = Country_Demand_projected_SSP_Incl_Losses.copy()

Node_Demand_SSP_projected_Incl_Losses = pd.merge(Node_Demand_2015[['Country' ,
                                                                   'PLEXOS_Countries' ,
                                                                   'PLEXOS_Nodes' , 
                                                                   'Share_%_Country_Demand']] , 
                                                 Node_Demand_SSP_projected_Incl_Losses  , 
                                                 left_on = 'Country' , 
                                                 right_index = True)

#Downscales projected country level demand by the 2015 shares of sub-country nodes as proxy.
for x in Years_List_5:
     Node_Demand_SSP_projected_Incl_Losses[x] = (Node_Demand_SSP_projected_Incl_Losses[x] * 
                                                 Node_Demand_SSP_projected_Incl_Losses['Share_%_Country_Demand']
                                                ).round(2).astype(float)
        
Node_Demand_SSP_projected_Incl_Losses = Node_Demand_SSP_projected_Incl_Losses.set_index('PLEXOS_Nodes')

#Filters data for the to-be modelled counties.
Node_Demand_SSP_projected_Incl_Losses = Node_Demand_SSP_projected_Incl_Losses[(Node_Demand_SSP_projected_Incl_Losses.
                                                                                 index.isin(Spatial_Mapping_Node.index))]

#Adds unit.
Node_Demand_SSP_projected_Incl_Losses.insert(loc = 3 , 
                                             column = 'Unit' , 
                                             value = 'GWh')

Node_Demand_SSP_projected_Incl_Losses[48:49]


# ### Interpolates 5-yearly values to yearly values and determines final electricity demand per node per year

# In[28]:


#Adds missing years to demand dataframe
for x in Years_List_Interpolation:
    
    Node_Demand_SSP_projected_Incl_Losses[x] = np.NaN

#Filters columns with numerical dtypes.
Node_Demand_SSP_projected_Incl_Losses_Years = Node_Demand_SSP_projected_Incl_Losses[Years_List]

#Linearly interpolates values between years.
Node_Demand_SSP_projected_Incl_Losses_Years = Node_Demand_SSP_projected_Incl_Losses_Years.interpolate(method = 'linear' , 
                                                                                                      axis = 1)

#Merges dataframes.
Node_Demand_SSP_projected_Incl_Losses = pd.merge(Node_Demand_SSP_projected_Incl_Losses[['Country' , 
                                                                                        'PLEXOS_Countries' , 
                                                                                        'Share_%_Country_Demand' , 
                                                                                        'Unit']] , 
                                                 Node_Demand_SSP_projected_Incl_Losses_Years , 
                                                 left_index = True , 
                                                 right_index = True)

#Saves the dataframe as csv.
Node_Demand_SSP_projected_Incl_Losses.to_csv(os.path.join(input_dir,'Final_Electricity_Demand_Nodal_Yearly.csv'))

Node_Demand_SSP_projected_Incl_Losses[48:49]
    
# ### Determines hourly peak demand per node per year

# In[29]:


#Calculates 2015 hourly peak demand per node in MWh
Node_Peak_Demand_2015 = Import_Hourly_Demand_2015.drop(columns = ['Datetime'])

Node_Peak_Demand_2015 = pd.DataFrame(Node_Peak_Demand_2015.max()).reset_index()

Node_Peak_Demand_2015.columns = ['PLEXOS_Nodes' , 
                                 'Node_Peak_Demand_2015']

#Calculates 2015 peak demand / total demand ratio
Node_Peak_Demand_2015 ['Ratio_Peak/Total_Demand_2015'] = Node_Peak_Demand_2015['Node_Peak_Demand_2015'
                                                                               ] / Node_Demand_2015['Node_Demand_2015']
#Merges dataframes.
Node_Peak_Demand_SSP_projected = pd.merge(Node_Peak_Demand_2015[['PLEXOS_Nodes' , 
                                                                 'Ratio_Peak/Total_Demand_2015']] , 
                                          Node_Demand_SSP_projected_Incl_Losses , 
                                          left_on = 'PLEXOS_Nodes' , 
                                          right_index = True , 
                                          how = 'right'
                                          ).set_index('PLEXOS_Nodes'
                                                     ).drop(columns = {'Share_%_Country_Demand'})

#Calculates projected hourly peak demand by using the relative 2015 peak demand as proxy. The peak to total demand ratio can be
#adjusted by changing the <Peak_Ratio_Factor>.
Peak_Ratio_Factor = 1

for x in Years_List:
    
     Node_Peak_Demand_SSP_projected[x] = (Node_Peak_Demand_SSP_projected[x] * 
                                          Peak_Ratio_Factor * 
                                          Node_Peak_Demand_SSP_projected['Ratio_Peak/Total_Demand_2015'] * 
                                          1000).round(2)

#Adds unit.
Node_Peak_Demand_SSP_projected['Unit'] = 'MW'

#Saves the dataframe as csv.
Node_Peak_Demand_SSP_projected.to_csv(os.path.join(input_dir,'Final_Electricity_Peak_Demand_Nodal_Yearly.csv'))
        
Node_Peak_Demand_SSP_projected.iloc[48:49]

model_start_year = parsed_yaml_file.get('startYear')
model_end_year = parsed_yaml_file.get('endYear')

#Format demand projections
with open(os.path.join(output_dir, 'SpecifiedAnnualDemand.csv'),'w') as f:
    f.write('REGION,FUEL,YEAR,VALUE\n')
    for x in Node_Demand_SSP_projected_Incl_Losses.index:
        if len(x) == 6:
            FUEL = 'ELC' + x[3:6] + 'XX02'
        if len(x) == 9:
            FUEL = 'ELC' + x[3:6] + x[7:9] + '02'
        for year in range(model_start_year,model_end_year+1):
            f.write('GLOBAL,' + str(FUEL) + ',' + str(year) + ',' + str(Node_Demand_SSP_projected_Incl_Losses.at[x, year]*(0.0036)) + '\n')
