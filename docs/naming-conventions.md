# Naming Conventions

This page describes in detial the naming conventios used. This is a helpful resource for when modifying the configuration file

## Power Technology Codes 

Each of the thirteen modelled power generation technologies are listed in the table below. Each power generation technology is given by a unique three letter code

| Technology                 | Code |
|----------------------------|------|
| Biomass                    | BIO  |
| Combined Cycle Natural Gas | CCG  |
| Coal                       | COA  |
| Cogeneration               | COG  |
| Concentrated Solar Power   | CSP  |
| Geothermal                 | GEO  |
| Hydroelectric              | HYD  |
| Open Cycle Natural Gas     | OCG  |
| Oil                        | OIL  |
| Other                      | OTH  |
| Petroleum                  | PET  |
| Solar Photovoltaic         | SPV  |
| Nuclear                    | URN  |
| Wave                       | WAV  |
| Waste                      | WAS  |
| Offshore Wind              | WOF  |
| Onshore Wind               | WON  |

## Commodity Codes 

Each of the power generation operate on a given commodity. All modelled commodites have a unique
three letter code, given in the table below. Note, that some technologies 
share commodities, such as Combined Cycle Natural Gas and Open Cycle Natural Gas 
both operating on GAS. 

| Commodity                  | Code |
|----------------------------|------|
| Biomass                    | BIO  |
| Natural Gas                | GAS  |
| Coal                       | COA  |
| Geothermal                 | GEO  |
| Hydroelectric              | HYD  |
| Oil                        | OIL  |
| Other                      | OTH  |
| Petroleum                  | PET  |
| Solar                      | SPV  |
| Nuclear                    | URN  |
| Wave                       | WAV  |
| Waste                      | WAS  |
| Offshore Wind              | WOF  |
| Onshore Wind               | WON  |

## Spatial Codes

### Country Codes 

OSeMOSYS Global can model a total of **163 countries**. All countries are described using their three 
letter codes, which can be found [here](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3). A few 
examples are given in the table below

| Country           | Code |
|-------------------|------|
| Netherlands       | NLD  |
| Great Britain     | GBR  |
| Canada            | CAN  |
| India             | IND  |

### Regional Node Codes

Some countries are further divided into regional nodes, for a total of **265 nodes**.
These nodes are represented by adding 2 characters to the country codes. If a country is 
represented by a single note, a XX placeholder is added. Some exaples are given in the table below.

x x x _ _ : First 3 characters (01-03) represent **country**

_ _ _ x x : Last 2 characters (04-05) represent **node**

| Country         | Node              | Code   |
|-----------------|-------------------|--------|
| United States   | California        | USACA  |
| Canada          | British Columbia  | CANBC  |
| Bangladash      | n/a               | BGDXX  |


### Fossil fuels and uranium:

x x x _ _ _ : First 3 characters (01-03) represent **fuel**

_ _ _ x x x : Last 3 characters (04-06) represent source; 
**country code** for **domestic** or '**INT**' for '**international**'. 

### Other fuels, including renewables: 8 characters 

^^^ _ _ _ _ _ : First 3  characters (01-03) represent **fuel**

_ _ _ ^^^ _ _ _ : Next 3 characters (04-06) represent **country code**

_ _ _ ^^^ _ _ _ : Last 2 characters (07-08) represent **node**

### Electricity:
### Electricity: 10 characters
^^^ _ _ _ _ _ _ _ : First 3 characters (01-03) are '**ELC**'

_ _ _ ^^^ _ _ _ _ : Next 3 characters (04-06) represent the **country code**

_ _ _ _ _ _ ^^ _ _ : Next 2 characters (07-08) represent the **node** 

_ _ _ _ _ _ _ _ ^^: Last 2 characters (09-10) represent the **level**; 
'**01**' for the secondary level (from powerplants to transmission) and '**02**' for the tertiary level (from transmission to distribution)

### Import and export

### Transmission and distribution

