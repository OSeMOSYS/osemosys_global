"""Holds lists of data files to be created."""

import pandas as pd 
from itertools import chain

POWER_PLANT_FILES = [
    'CapitalCost.csv',
    'FixedCost.csv',
    'CapacityToActivityUnit.csv',
    'OperationalLife.csv',
    'TotalAnnualMaxCapacityInvestment.csv',
    'TotalAnnualMinCapacityInvestment.csv',
    'TotalTechnologyModelPeriodActivityUpperLimit.csv',
    'FUEL.csv',
    'InputActivityRatio.csv',
    'OutputActivityRatio.csv',
    'MODE_OF_OPERATION.csv',
    'REGION.csv',
    'ResidualCapacity.csv',
    'TECHNOLOGY.csv',
    'YEAR.csv'
    ]

TIMESLICE_FILES = [
    'CapacityFactor.csv',
    'TIMESLICE.csv',
    'SpecifiedDemandProfile.csv',
    'YearSplit.csv'
    ]

VARIABLE_COST_FILES = [
    'VariableCost.csv'
    ]

DEMAND_FILES = [
    'SpecifiedAnnualDemand.csv'
    ]

EMISSION_FILES = [
    'EmissionActivityRatio.csv',
    'EmissionsPenalty.csv',
    'EMISSION.csv'
]

MAX_CAPACITY_FILES = [
    'TotalAnnualMaxCapacity.csv'
]

USER_CAPACITY_FILES = [
    'TotalAnnualMinCapacityInvestment.csv',
    'TotalAnnualMaxCapacityInvestment.csv'
]

DEMAND_FIGURES = [
    'South America',
    'Oceania',
    'North America',
    'Europe',
    'Asia',
    'Africa'
]

# Files to be created by file_check rule
INPUT_FILES = pd.read_csv('./resources/otoole_files.csv')['inputs'].to_list()
CHECK_FILES = [f'{f}.csv' for f in INPUT_FILES]
GENERATED_FILES = list(chain(
    POWER_PLANT_FILES, 
    TIMESLICE_FILES, 
    VARIABLE_COST_FILES, 
    DEMAND_FILES, 
    EMISSION_FILES, 
    MAX_CAPACITY_FILES))
for csv in GENERATED_FILES:
    CHECK_FILES.remove(csv) # remove .csv
CHECK_FILES.remove('default_values.csv')

