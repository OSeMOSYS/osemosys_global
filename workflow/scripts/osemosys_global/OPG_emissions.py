"""Generates emission set, emission activity ratio, and emission penalty"""

import logging
import pandas as pd
from pathlib import Path
from OPG_configuration import ConfigFile, ConfigPaths


# Logging formatting 
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

# Constant for tech to fuel emission type mapping 
_TECH_TO_FUEL = {
    'BIO':'Bagasse',
    'COA':'Lignite Coal',
    'COG':'Lignite Coal',
    'OCG':'Natural Gas',
    'CCG':'Natural Gas',
    'PET':'Motor Gasoline',
    'OIL':'Crude Oil',
    'OTH':'Natural Gas'
}

# Emission name 
_EMISSION = 'CO2'

def main():
    """Assigns CO2 equivalent values to each technology over all its modes 
    of operation. For technologies that do not have emissions, it assigns a 
    value of zero. A global emission penalty value is applied and the 
    emission type (CO2) is written to the emission set. 
    """

    # CONFIGURATION PARAMETERS

    config_paths = ConfigPaths()
    config = ConfigFile('config')  
    
    output_data_dir = config_paths.output_data_dir
    emission_penalty = config.get('emission_penalty') # M$/MT

    # ASSIGN EMISSION ACTIVITY RATIOS

    df_ear = get_ear(_EMISSION)
    df_ear.to_csv(Path(output_data_dir, 'EmissionActivityRatio.csv'), index=False)
    logging.info('Successfully generated emission activity ratio')

    # ASSIGN EMISSION PENALTY 

    df_emission_penalty = get_emission_penalty(_EMISSION, emission_penalty)
    df_emission_penalty.to_csv(Path(output_data_dir, 'EmissionsPenalty.csv'), index=False)
    logging.info('Successfully generated emission penalty')

    # ASSIGN EMISSION 

    df_emission = pd.DataFrame([_EMISSION], columns=['VALUE'])
    df_emission.to_csv(Path(output_data_dir, 'EMISSION.csv'), index=False)
    logging.info('Successfully generated emission set')


def get_co2_emission_factors():
    """Gets co2 emission factors for diferent fuels. 

    Reads in a file containing co2, ch4, n2o emission factors and global 
    warming potentials for various fuel types. This function performs unit 
    conversions to convert everyting to MegaTonnes per PetaJoule and collapses
    the different emission types to a single co2 equivalent value for each 
    fuel. 

    Returns:
        Dictionary holding fuel type as key and co2 factor as value in units 
        of MT/PJ 

    Example: 
        co2_factors = get_co2_emission_factors()
        co2_factors['Natural Gas'] 
        -> 0.0503
    """

    # Source for all emission factors comes from: 
    # https://www.epa.gov/sites/default/files/2018-03/documents/emission-factors_mar_2018_0.pdf

    # Configuration parameters
    config_paths = ConfigPaths()

    # Read in emission factors 
    input_data_dir = config_paths.input_data_dir
    df_raw = pd.read_csv(Path(input_data_dir,'emission_factors.csv'))
    df_raw = df_raw.drop([0]).reset_index(drop=True) # drop units row

    # Convert co2 factors from kg/mmbtu to MT/PJ 
    # kg/mmbtu * 1mmbtu/1.05GJ * 1000000GJ / PJ * 1T/1000kg * 1MT/1000000T
    # Multiply by global warming potential to get co2_eq
    co2 = df_raw['co2_factor'].astype(float) * (1/1055) * df_raw['co2_gwp'].astype(float)

    # Convert ch4 and n2o factors from g/mmbtu to MT/PJ 
    # kg/mmbtu * 1mmbtu/1.05GJ * 1000000GJ / PJ * 1T/1000000g * 1MT/1000000T
    # Multiply by global warming potential to get co2_eq
    ch4 = df_raw['ch4_factor'].astype(float) * (1/1055000) * df_raw['ch4_gwp'].astype(float)
    n2o = df_raw['n2o_factor'].astype(float) * (1/1055000) * df_raw['n2o_gwp'].astype(float)

    # Find total CO2 equivalent
    data = {'co2':co2, 'ch4':ch4, 'n2o':n2o} 
    df = pd.DataFrame(data).set_axis(df_raw['FUEL TYPE']) 
    df['co2_eq'] = round(df.sum(axis=1),4).astype(float)

    return df.set_index(df.index).to_dict()['co2_eq']

def get_ear(emission):
    """Creates emission activity ratio dataframe.
    
    This function reads in an existing input activity ratio parameter file
    and removes the fuel and year columns. This leaves a dataframe with info
    on when all technologies are allowed to operate over the model horizon. 
    A column is added in to hold the emission type and emission activity ratio
    based. 

    Args: 
        emission: string describing the emission type (ie. 'CO2')

    Returns:
        df: Dataframe describing emission activity ratio. Dataframe headers 
            shown below 

            REGION | TECHNOLOGY | EMISSION | MODE_OF_OPERATION | YEAR | VALUE
            ----------------------------------------------------------|-----
               x   |      x     |    x     |        x          |  x   |  x
    """

    # CONFIGURATION PARAMETERS
    config_paths = ConfigPaths()
    output_data_dir = config_paths.output_data_dir

    # GET EMISSION FACTORS 

    co2_factors = get_co2_emission_factors()

    # GET INFO FROM INPUT ACTIVITY RATIO

    df_iar = pd.read_csv(Path(output_data_dir, 'InputActivityRatio.csv'))
    df = df_iar.drop(['FUEL', 'VALUE'], axis=1)

    # ADD IN EMISSION COLUMN

    df['EMISSION'] = emission
    
    # ADD MAPPING OF TECHNOLOGY TO EMISSION ACTIVITY RATIO

    techs = pd.Series(df['TECHNOLOGY'].str[3:6])
    fules = techs.map(_TECH_TO_FUEL)
    df['VALUE'] = fules.map(co2_factors)
    df['VALUE'] = df['VALUE'].fillna(0)
    df = df[[
        'REGION', 
        'TECHNOLOGY', 
        'EMISSION', 
        'MODE_OF_OPERATION', 
        'YEAR', 
        'VALUE']]
    
    return df

def get_emission_penalty(emission, penalty): 
    """Creates emission penalty dataframe. 

    The emission penalty is applied at a global geographical level. All regions
    and subregions have the same penalty. 

    Args: 
        emission: string describing the emission type (ie. 'CO2')
        penalty: emission penalty in M$/MT

    Returns: 
        df: Dataframe describing emission penalty. Dataframe headers shown 
            below

            REGION | EMISSION | YEAR | VALUE
            ---------------------------------
               x   |    x     |  x   |   x
    """

    # CONFIGURATION PARAMETERS
    config = ConfigFile('config')
    start_year = config.get('startYear')
    end_year = config.get('endYear')
    region = config.get('region')

    # GENERATE DATA
    
    data = []
    for year in range(start_year, end_year+1):
        data.append([region, emission, year, penalty])
    df = pd.DataFrame(data, columns=['REGION','EMISSION','YEAR','VALUE'])
    return df

if __name__ == '__main__':
    main()
