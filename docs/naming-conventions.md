# Naming Conventions

This page describes the detials of the naming conventions used in OSeMOSYS 
Global. This is a helpful resource for when modifying the configuration file 
or interpreting the results. 

## Reference Energy System

The schematic below shows a high level overview of the reference energy system 
(RES). A full RES example for a single node can be found [here](#full-res-example)

:::{seealso}
OSeMOSYS' 
[documentation](https://osemosys.readthedocs.io/en/latest/manual/Create%20a%20model%20in%20OSeMOSYS.html#mapping-the-res-of-atlantis) 
on Reference Energy Systems
:::

## OSeMOSYS Global Structure

OSeMOSYS models are built through connecting commodities (ie. fules) and 
energy conversion technologies. This page will walk through how OSeMOSYS 
Global connects and names each technology and commodity. 

### Spatial Codes

#### Country Codes 

OSeMOSYS Global can model a total of **163 countries**. All countries are 
described using their three letter codes, which can be found 
[here](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3). A few 
examples are given in the table below

| Country           | Code |
|-------------------|------|
| Netherlands       | NLD  |
| Great Britain     | GBR  |
| Canada            | CAN  |
| India             | IND  |

#### Regional Node Codes

Some countries are further divided into regional nodes, for a total of **265 
nodes**. These nodes are represented by adding 2 characters to the country 
codes. If a country is represented by a single note, a XX placeholder is added. 

| Character Example | Location | Length | Description        |
|-------------------|----------|--------|--------------------|
| X X X _ _         | (01-03)  | 5      | Country Code       |
| _ _ _ X X         | (04-05)  | 5      | Regional Node Code |

Examples of spatial codes are given below

| Country         | Node              | Code   |
|-----------------|-------------------|--------|
| United States   | California        | USACA  |
| Canada          | British Columbia  | CANBC  |
| Bangladash      | n/a               | BGDXX  |

### Acronyms 

Each unique technology and commodity in OSeMOSYS Global has a three-character
code used to identify it. These three letter codes are not directly modelled, 
instead they are used to build the full technologies and commodity codes used
in OSeMOSYS Global. 

#### Technology Acronyms

Thirteen power generation technologies are modelled in OSeMOSYS Global. The 
acronyms used to identify each technology are given below. 

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

#### Commodity Acronyms

Each of the power generation technologies operate on a given commodity. The 
acronyms used to identify these commodities are given in the table below. Note, 
that some technologies share commodities, such as Combined Cycle Natural Gas 
and Open Cycle Natural Gas both operating on the commodity GAS. 

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

### Technology Codes

Technologies are responsible for converting one energy carrier into another. 
Each technology can have unique costs, efficiencies, and capacity/generation 
limits. There are four types of technologies in OSeMOSYS Global; mining 
technologies, power generation technoloiges, trasmission technologies, and 
trade technologies. 

#### Mining Technology Codes

Mining technologies are responsible for introducing raw commodities into the 
model. These technologies are either 9 or 11 characters long, depending on 
the type of resource the technology mines. Technologies that introduce 
non-tradable renewable resources into the model are 11 characters long, while 
technologies that introduce tradable commodities are 9 letters long. The
difference comes from renewable resources being tracked at a **nodal** level, 
while physical resources being tracked at a **country** level.

| Character Example         | Location | Length | Description |
|---------------------------|----------|--------|-------------|
| M I N _ _ _ _ _ _ _ _     | (01-03)  | 11     | Mining tradable commodity |
| R N W _ _ _ _ _ _         | (01-03)  | 9      | Mining non-tradable commodity |
| _ _ _ X X X _ _ _ ( _ _ ) | (04-06)  | 9 / 11 | Technology code |
| _ _ _ _ _ _ X X X ( _ _ ) | (07-09)  | 9 / 11 | Country code or `INT` for International* |
| R N W _ _ _ _ _ _ X X     | (10-11)  | 11     | Regional node code |

(*) See trading technologies for more information

Examples of mining technologies are given in the table below

| Code        | Description |
|-------------|------|
| MINGASCAN   |   |
| MINGASINT   |   |
| RNWWNSCANBC |   |

#### Power Generation Technology Codes

Power generation technologies are responsible for converting raw commodities
into electricity transmission commodities. These technologies are 13 
characters long and include spatial and technology identifiers. The makeup of 
these codes are given in the table below. 

| Character Example         | Location | Length | Description |
|---------------------------|----------|--------|-------------|
| P W R _ _ _ _ _ _ _ _ _ _ | (01-03)  | 13     | Power generator |
| _ _ _ X X X _ _ _ _ _ _ _ | (04-06)  | 13     | Technology code |
| _ _ _ _ _ _ X X X _ _ _ _ | (07-09)  | 13     | Country code |
| _ _ _ _ _ _ _ _ _ X X _ _ | (10-11)  | 13     | Regional node code |
| _ _ _ _ _ _ _ _ _ _ _ 0 0 | (12-13)  | 13     | Technology that can **not** be invested in |
| _ _ _ _ _ _ _ _ _ _ _ 0 1 | (12-13)  | 13     | Technology that **can** be invested in |

All power generation technologies will at minimum be represented with a `01` 
at the end of the code, signifying the technology can be invested in (if the 
total max capacity limits allow it). Only some technologies will also be 
represented by coes that end with a `00`. These will represent historical 
technologies that exist at the start of the model but can no longer be invested
in, such as thermal generators built a number of years ago with poor 
efficiencies. 

Examples of power generation technologies are given in the table below

| Code          | Description |
|---------------|-------------|
| PWRURNCANBC01 |   |
| PWROCGINDNE00 |   |
| PWRSPVUSACA01 |   |

#### Transmission Technology Codes

Transmission technologies are responsible for coverting a regional node's 
transmission electricity into end use electricity, and recieving traded 
electricity. It acts as a dummy aggregation technology. The makeup of this 
code is given in the table below. 

| Character Example     | Location | Length | Description |
|-----------------------|----------|--------|-------------|
| P W R _ _ _ _ _ _ _ _ | (01-03)  | 11     | Power generator |
| _ _ _ T R N _ _ _ _ _ | (04-06)  | 11     | Transmission technology |
| _ _ _ _ _ _ X X X _ _ | (07-09)  | 11     | Country code |
| _ _ _ _ _ _ _ _ _ X X | (10-11)  | 11     | Regional node code |

Example transmisison technologies are given below. 

| Code        | Description |
|-------------|-------------|
| PWRTRNCANBC |   |
| PWRTRNINDNE |   |
| PWRTRNUSACA |   |

### Trading Technology Codes

:::{seealso}
Out description of how resource trading functions in OSeMOSYS Global
:::

Trading technologies are responsible for trading electricity. While physical 
commodities can also be traded, these are handeled through trading internationl
fules. Trading technologies connect from one nodes end use electricity to another 
nodes transmission electricity. The table below summarises the codes. 

| Character Example         | Location | Length | Description |
|---------------------------|----------|--------|-------------|
| T R N _ _ _ _ _ _ _ _ _ _ | (01-03)  | 13     | Transmission technology |
| _ _ _ X X X _ _ _ _ _ _ _ | (04-06)  | 13     | Country to trade **from** |
| _ _ _ _ _ _ X X _ _ _ _ _ | (07-09)  | 13     | Regional node to trade **from** |
| _ _ _ _ _ _ _ _ X X X _ _ | (10-11)  | 13     | Country to trade **to** |
| _ _ _ _ _ _ _ _ _ _ _ X X | (12-13)  | 13     | Regional node to trade **to** |

### Commodity Codes 

Commodities (or fuel) carry energy and flow into and out of technologies in OSeMOSYS. 
There are three classifications of commodities in OSeMOSYS Global, raw fuel, 
transmission level fuel, and end use fuel. 

#### Raw Fuel Codes

Raw fules in OSeMOSYS Global can be three, six, or nine characters long.
Depending on the fuels purpose, the country and/or region may or may not need 
to be tracked. 


<!-- 
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

### Transmission and distribution -->

## Full Example

The schematic below shows the reference energy system for the node of 
British Columbia, Canada. All nodes follow the similar structure. 

![RES](_static/res.jpg "RES")