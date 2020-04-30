# Naming convention
## Country codes: 3 characters
Total of **166 countries**.

Examples: 
**NLD** for the Netherlands, **GBR** for Great Britain, 
**CAN** for Canada, and **IND** for India

#### Nodes: 
Some countries are further divided into regional **nodes**. Total of **265 nodes**.
These nodes are represented by adding 2 characters to the country codes

Examples: 
**USACA** for California, USA, **CANBC** for British Columbia, Canada

Note:
'**XX**' is added as a placeholder for countries which are not sub-divided into nodes. 

## Commodity codes: 6/8 characters
### Fossil fuels and uranium: 6 characters
^^^ _ _ _ : First 3 characters (01-03) represent **fuel**

_ _ _ ^^^ : Last 3 characters (04-06) represent source; 
**country code** for **domestic** or '**INT**' for '**international**'. 

### Other fuels, including renewables: 8 characters 
^^^ _ _ _ _ _ : First 3  characters (01-03) represent **fuel**

_ _ _ ^^^ _ _ _ : Next 3 characters (04-06) represent **country code**

_ _ _ ^^^ _ _ _ : Last 2 characters (07-08) represent **node**

### Electricity:
=======
### Electricity: 10 characters
^^^ _ _ _ _ _ _ _ : First 3 characters (01-03) are '**ELC**'

_ _ _ ^^^ _ _ _ _ : Next 3 characters (04-06) represent the **country code**

_ _ _ _ _ _ ^^ _ _ : Next 2 characters (07-08) represent the **node** 

_ _ _ _ _ _ _ _ ^^: Last 2 characters (09-10) represent the **level**; 
'**01**' for the secondary level (from powerplants to transmission) and '**02**' for the tertiary level (from transmission to distribution)

## Technology codes

The following technologies exist in the model:
**COA** represents coal
**COG** represents coal gas
**CCG** represents natural gas combined cycle
**OCG** represents natural gas open cycle
**PET** represents petroleum
**URN** represents uranium
**OIL** represents oil and oil products
**OTH** represents other fuels

The following renewable fuels are included:
**BIO** represents biomass
**GEO** represents geothermal
**HYD** represents hydro
**SPV** represents solar
**CSP** represents solar
**WAS** represents waste recovery
**WAV** represents wave
**WON** represents onshore wind
**WOF** represents offshore wind

### Powerplants: 13 characters

### Fuel codes

The following fuels are included:
**COA** represents coal
**COG** represents coal gas
**GAS** represents natural gas
**PET** represents petroleum
**URN** represents uranium
**OIL** represents oil and oil products
**OTH** represents other fuels

The following renewable fuels are included:
**BIO** represents biomass
**GEO** represents geothermal
**HYD** represents hydro
**SOL** represents solar
**WAS** represents waste recovery
**WAV** represents wave
**WIN** represents wind

### Import and export

### Transmission and distribution

