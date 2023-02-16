# Indonesia

The first application of FEO is the development of an electricity systems model 
for Indonesia. 

## Model scope

### Spatial resolution
The model represents 34 provinces of Indonesia across 7 regions, as shown in the map below:
 - Sumatra [10 provinces]
 - Java [6 provinces]
 - Kalimantan [5 provinces]
 - Nusa Tenggara [3 provinces]
 - Sulawesi [6 provinces]
 - Maluku [2 provinces]
 - Papua [2 provinces]

![IDN-Provinces](figures/Indonesia_provinces_english.png "IDN-Provinces")

### Temporal resolution
Each year is divided into 6 'Seasons' [S1-S6]: 
Each 'Season' is further divided into 12 'Daily Time Brackets' 
Together, there are 72 representative 'timeslices' in the model. The temporal resolution is the same for the entire model period. 


### Model horizon
Base year - 2021 \
End year - 2050

## Key assumptions

### Discount rates
The model includes two types of discount rates (DR): 'social' and 'financial'. 
The social DR is applied across the entire model and represents the relative 
weighting of present and future costs and benefits. A low social DR weights the present and the future more similarly than a high DR. The financial DR is technology-specific and represents the weighted average cost of capital (WACC) for a given technology (e.g. power plant). The model assumes a value of **10%** for both the social and financial discount rates. The latter is based on the [IEA's Cost of Capital Observatory](https://www.iea.org/data-and-statistics/data-tools/cost-of-capital-observatory)

### Reserve Margin

- Reserve margin 60% -> 35%
- 

## Data
### Technology costs

- Capital costs (USD<sub>2020</sub>/kW, except for *Battery storage* which is in USD <sub>2020</sub>/kWh)

| Technology          | 2020    | 2030    | 2040    | 2050    | Source  |
|---------------------|---------|---------|---------|---------|---------|
| Battery storage	  |         |         |         |         | [NREL](https://www.nrel.gov/docs/fy21osti/79236.pdf, 'NREL')  |
| Coal                |         |         |         |         |         |
| Gas - CCGT            |         |         |         |         |         |
| Gas - OCGT            |         |         |         |         |         |
| Hydropower          |         |         |         |         |         |
| Biomass             |         |         |         |         |         |
| Solar photovoltaic |         |         |         |         |         |
| Wind - Onshore          |         |         |         |         |         |
| Wind- Offshore          |         |         |         |         |         |
| CSP          |         |         |         |         |         |
| Geothermal          |         |         |         |         |         |
| Diesel          |         |         |         |         |         |


### Renewable Energy Profiles

- renewables.ninja

### Renewable Energy Potentials

- GEO: Volcanostratigraphy of Batukuwung-Parakasak Geothermal
Area, Serang Regency, West Java. 
- SPV, WON, HYD: Beyond 443 GW

| Renewable source  | Potential | Source    |
|-------------------|-----------|-----------|
| Solar PV          |           |           |
| Wind - Onshore    |           |           |
| Hydropower        |           |           |
| Geothermal        |           |           |
| Biomass           |           |           |

### Energy demand projections

- Own calculations

### Fuel Prices


| Fuel          |Unit       | 2020    | 2030    | 2040    | 2050    | Source  |
|---------------|-----------|---------|---------|---------|---------|---------|
| Coal	        | USD/mt    |         |         |         |         |         |
| Natural gas   | USD/mmbtu |         |         |         |         |         |
| Oil           | USD/bbl   |         |         |         |         |         |

### Electricity interconnectors


## Scenarios

The model was used to explore three scenarios: *Current Policies [CP]*, 
*Least-cost [LC]*, and *Net-Zero [NZ]*. The scenarios represent alternate 
pathways for the expansion of Indonesia's electricity system. Each scenario 
consists of a set of assumptions and constraints, as detailed below:

### Current policies

This scenario includes all implemented policies related to the expansion of 
Indonesia's electricity system as well as power plants under construction. 
The policies included are: 

And the future power plants included are:

### Least-cost


### Net-zero


## Results


## Planned improvements

- Interconnector expansion plans
- Fossil fuel price projections
- Plant-specific efficiencies
- Hydropower capacity factor by plant / node
- Technology-specific discount rates

