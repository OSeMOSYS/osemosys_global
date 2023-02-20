Indonesia
=========

The first application of FEO Global is the development of an openly available 
electricity systems model for Indonesia. This model is used to explore 
transition pathways to a net-zero electricity system. As with all energy system 
models, the inputs include a range of datasets and assumptions (e.g. technolgy 
cost projections, discount rates). All of these inputs are described here in 
order to allow for the model to be reviewed, re-run, and re-purposed.

Model scope
-----------
The model aims to represent the electricity system of Indonesia as accurately as
possible, subject to constraints on data and computation time. The main aspects 
that improve the accuracy of the model's repersentation of Indonesia's 
electricity system are the spatial and temporal resolution.

Spatial resolution
..................

The model represents 34 provinces of Indonesia across 7 regions, 
as shown in the table and map below:

.. csv-table:: 
   :file: tables/provinces.csv
   :widths: 50, 50, 50, 50
   :header-rows: 1
   :align: center

.. image:: figures/Indonesia_provinces_english.png


Temporal resolution
...................
Each year is divided into 6 'Seasons' [S1-S6]: 

.. csv-table:: 
   :file: tables/seasons_summary.csv
   :widths: 50, 50
   :header-rows: 1
   :align: center

Each 'Season' is further divided into 12 'Daily Time Brackets':

.. csv-table:: 
   :file: tables/dailytimebrackets_summary.csv
   :widths: 50, 50
   :header-rows: 1
   :align: center

Together, there are 72 representative 'timeslices' in the model. The temporal 
resolution is the same for the entire model period. 


Model horizon
.............
Base year - 2021

End year - 2050

Key assumptions
---------------

Discount rates
..............
The model includes two types of discount rates (DR): 'social' and 'financial'. 
The social DR is applied across the entire model and represents the relative 
weighting of present and future costs and benefits. A low social DR weights the 
present and the future more similarly than a high DR. The financial DR is 
technology-specific and represents the weighted average cost of capital (WACC) 
for a given technology (e.g. power plant). The model assumes a value of 
**10%** for both the social and financial discount rates. The latter is based 
on the `IEA Cost of Capital Observatory <iea_wacc_>`_ 


.. _iea_wacc: https://www.iea.org/data-and-statistics/data-tools/cost-of-capital-observatory

Reserve Margin
..............

- Reserve margin 60% -> 35%
- 

Data
----

Technology costs
................

`NREL <https://www.nrel.gov/docs/fy21osti/79236.pdf>`_ 

.. csv-table:: Technology cost projections (Capital)
   :file: tables/technology_costs_capital.csv
   :widths: 75, 50, 50, 50, 50, 50
   :header-rows: 1
   :align: center

.. csv-table:: Technology cost projections (Fixed)
   :file: tables/technology_costs_fixed.csv
   :widths: 75, 50, 50, 50, 50, 50
   :header-rows: 1
   :align: center


Renewable Energy Profiles
.........................
- renewables.ninja

Renewable Energy Potentials
...........................
- GEO: Volcanostratigraphy of Batukuwung-Parakasak Geothermal \
  Area, Serang Regency, West Java. 
- SPV, WON, HYD: Beyond 443 GW
- HYD: https://www.hydropower.org/blog/indonesia-promotes-hydropower-to-create-the-demand-for-industrial-development#:~:text=The%20biggest%20hydropower%20potential%20is,Tenggara%2DMaluku%20is%201.1%20GW.


.. csv-table:: Renewable energy potential
   :file: tables/re_potentials_summary.csv
   :widths: 75, 50, 50, 50, 50, 50
   :header-rows: 1
   :align: center

Energy demand projections
.........................
- Own calculations

Fuel Prices
...........

.. csv-table:: Fuel price projections
   :file: tables/fuel_prices.csv
   :widths: 75, 75, 50, 50, 50, 50
   :header-rows: 1

Electricity interconnectors
...........................


Scenarios
---------

The model was used to explore three scenarios: *Current Policies [CP]*, 
*Least-cost [LC]*, and *Net-Zero [NZ]*. The scenarios represent alternate 
pathways for the expansion of Indonesia's electricity system. Each scenario 
consists of a set of assumptions and constraints, as detailed below:

Current policies
................

This scenario includes all implemented policies related to the expansion of 
Indonesia's electricity system as well as power plants under construction. 
The policies included are: 

And the future power plants included are:


Least-cost
..........


Net-zero
........



Results
-------

Capacity expansion
..................

.. raw:: html
   :file: figures/TotalCapacityAnnual_BAU.html

Annual electricity generation mix
.................................

.. raw:: html
   :file: figures/GenerationAnnual_BAU.html

Hourly electricity generation mix
.................................

.. raw:: html
   :file: figures/GenerationHourly_BAU.html


Planned improvements
--------------------

* Interconnector expansion plans
* Fossil fuel price projections
* Plant-specific efficiencies
* Hydropower capacity factor by plant / node
* Technology-specific discount rates

Model code, data, and workflow
------------------------------

The entire workflow of FEO Global is available under an open license 
at `transition-zero/feo-esmod-osemosys`. 
In addition, it uses only publicly available data and open source solver (CBC). 


References
----------

* IEA 
* IRENA
* RUPTL 2021-2030
* Beyond 443 GW
* LTS-LCCP (Indonesia submission to UNFCCC)