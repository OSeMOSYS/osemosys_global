import pandas as pd
import os

from typing import Optional, Any

import logging

logger = logging.getLogger(__name__)

from read import(
    import_storage_build_rates,
    import_op_life,
    import_GESDB_project_data,
    import_GESDB_regional_mapping,
    import_iar_base,
    import_oar_base,
    import_capact_base,
    import_cap_cost_base,
    import_fix_cost_base,
    import_var_cost_base,
    import_op_life_base,
    import_max_cap_invest_base,
    import_min_cap_invest_base,
    import_res_cap_base,
    import_set_base
)

from constants import(
    GESDB_TECH_MAP,
    DURATION_TYPE,
    BUILD_YEAR,
    RETIREMENT_YEAR,
    INACTIVE_VARS,
    ACTIVE_VARS,
    NEW_VARS
    )

from activity import(
    activity_storage,
    create_storage_capacity_activity
    )

from costs import(
    set_storage_capex_costs,
    set_storage_operating_costs
    )
    
from operational_life import set_op_life_storage

from investment_constraints import cap_investment_constraints_sto

from user_defined_capacity import set_user_defined_capacity_sto

from residual_capacity import res_capacity_storage

from sets import(set_unique_storage_technologies, 
                 set_unique_technologies)

from technology_to_from_storage import (set_technology_to_storage,
                                        set_technology_from_storage
                                        )

from storage_level import set_storage_level_start

def main(
    build_rates: pd.DataFrame,
    unique_sto_techs: list,
    default_op_life: dict[str, int],
    gesdb_data: pd.DataFrame,
    gesdb_mapping: pd.DataFrame,
    iar_base: pd.DataFrame,
    oar_base: pd.DataFrame,
    capact_base: pd.DataFrame,
    cap_cost_base: pd.DataFrame,
    fix_cost_base: pd.DataFrame,
    var_cost_base: pd.DataFrame,    
    op_life_base: pd.DataFrame,
    max_cap_invest_base: pd.DataFrame,
    min_cap_invest_base: pd.DataFrame,
    res_cap_base: pd.DataFrame,
    tech_set_base: pd.DataFrame,
    fuel_set_base: pd.DataFrame,
    tech_capacity_storage: Optional[dict[str, list[Any]]] = None
):
    
    if not unique_sto_techs:
        
        logger.warning("No storage added to the system. Populate 'storage_parameters' in config file.")
        
        oar_storage = oar_base.copy()
        iar_storage = iar_base.copy()
        cap_activity_storage = capact_base.copy()
        cap_cost_storage = pd.DataFrame(columns=["REGION", "STORAGE", "YEAR", "VALUE"])
        cap_cost = cap_cost_base.copy()
        fix_cost_storage = fix_cost_base.copy() 
        var_cost_storage = var_cost_base.copy()
        op_life = op_life_base.copy()
        op_life_storage = pd.DataFrame(columns=["REGION", "STORAGE", "VALUE"])
        max_cap_invest_storage = max_cap_invest_base.copy()
        tech_set = tech_set_base.copy()
        storage_set = pd.DataFrame(columns=["VALUE"])
        tech_to_storage = pd.DataFrame(columns=["REGION","TECHNOLOGY","STORAGE","MODE_OF_OPERATION", "VALUE"])
        tech_from_storage = pd.DataFrame(columns=["REGION","TECHNOLOGY","STORAGE","MODE_OF_OPERATION", "VALUE"])
        res_cap = res_cap_base.copy()
        res_cap_storage = pd.DataFrame(columns=["REGION", "STORAGE", "YEAR", "VALUE"])
        storage_level_start = pd.DataFrame(columns=["REGION", "STORAGE", "VALUE"])
        tech_capacity_storage = None

    else:
    
        # Set storage set
        storage_set = set_unique_storage_technologies(fuel_set_base, unique_sto_techs)

        # Create TechnologyToStorage and TechnologyFromStorage
        tech_to_storage = set_technology_to_storage(storage_set, region_name)
        tech_from_storage = set_technology_from_storage(storage_set, region_name)    

        
        # Set capital costs for storage.
        cap_cost, cap_cost_storage = set_storage_capex_costs(storage_set, 
                                                            storage_parameters,
                                                            cap_cost_base,
                                                            start_year,
                                                            end_year,
                                                            region_name)
        
        # Set fixed and variable operating costs for storage.
        fix_cost_storage, var_cost_storage = set_storage_operating_costs(storage_set,
                                                                        storage_parameters,
                                                                        fix_cost_base,
                                                                        var_cost_base,
                                                                        start_year,
                                                                        end_year,
                                                                        region_name)
        
        # Set activity ratios for storage.
        iar_storage, oar_storage = activity_storage(storage_set, iar_base, oar_base, 
                                                    storage_parameters, start_year, 
                                                    end_year, region_name)
        
        # assign capacity to activity unit to storage.
        cap_activity_storage = create_storage_capacity_activity(storage_set, capact_base)
        
        # Set operational life for storage.
        op_life, op_life_storage = set_op_life_storage(storage_set, default_op_life, 
                                                    op_life_base, region_name)
        
        # Set annual capacity investment constraints.
        max_cap_invest_storage = cap_investment_constraints_sto(storage_set, 
                                                                max_cap_invest_base,
                                                                build_rates,
                                                                no_investment_techs, 
                                                                start_year, 
                                                                end_year, 
                                                                region_name)
        
        # Set residual capacity ('PWR') and residual capacity storage.
        res_cap, res_cap_storage = res_capacity_storage(gesdb_data, gesdb_mapping, 
                                                        res_cap_base, storage_existing, 
                                                        storage_planned, op_life_dict, 
                                                        storage_parameters,
                                                        GESDB_TECH_MAP, DURATION_TYPE,
                                                        BUILD_YEAR, RETIREMENT_YEAR,
                                                        INACTIVE_VARS, ACTIVE_VARS,
                                                        NEW_VARS, start_year,
                                                        end_year, region_name)
        
        storage_level_start = set_storage_level_start(storage_set, region_name)

        # Alter output csv's based on user defined capacities following user config.
        if tech_capacity_storage is not None:
            (max_cap_invest_storage, 
            min_cap_invest_storage, 
            res_cap,
            res_cap_storage,
            oar_storage,
            cap_cost,
            cap_cost_storage,
            fix_cost_storage, 
            var_cost_storage
            ) = set_user_defined_capacity_sto(
                tech_capacity_sto, 
                storage_parameters,
                default_op_life, 
                min_cap_invest_base, 
                max_cap_invest_storage, 
                res_cap,
                res_cap_storage,
                oar_storage,
                cap_cost,
                cap_cost_storage,
                fix_cost_storage, 
                var_cost_storage,
                start_year,
                end_year,
                region_name
                )  

        # get new additions to technology sets
        new_techs = set_unique_technologies(storage_set)
        tech_set = pd.concat([tech_set_base, new_techs])
        
    # OUTPUT CSV's
    
    oar_storage.to_csv(os.path.join(output_data_dir, "OutputActivityRatio.csv"), index=None)
    
    iar_storage.to_csv(os.path.join(output_data_dir, "InputActivityRatio.csv"), index=None)

    cap_activity_storage.to_csv(os.path.join(output_data_dir, "CapacityToActivityUnit.csv"), index = None)
    
    cap_cost_storage.to_csv(os.path.join(output_data_dir, "CapitalCostStorage.csv"), index = None)
    cap_cost.to_csv(os.path.join(output_data_dir, "CapitalCost.csv"), index = None)
    
    fix_cost_storage.to_csv(os.path.join(output_data_dir, "FixedCost.csv"), index = None)
    var_cost_storage.to_csv(os.path.join(output_data_dir, "VariableCost.csv"), index = None)
    
    op_life_storage.to_csv(os.path.join(output_data_dir, "OperationalLifeStorage.csv"), index = None)
    op_life.to_csv(os.path.join(output_data_dir, "OperationalLife.csv"), index = None)
    
    max_cap_invest_storage.to_csv(os.path.join(output_data_dir, 
                                            'TotalAnnualMaxCapacityInvestment.csv'),
                                        index = None)

    tech_set.to_csv(os.path.join(output_data_dir, "TECHNOLOGY.csv"), index = None)
    storage_set.to_csv(os.path.join(output_data_dir, "STORAGE.csv"), index = None)
    
    tech_to_storage.to_csv(os.path.join(output_data_dir, "TechnologyToStorage.csv"), index=None)
    tech_from_storage.to_csv(os.path.join(output_data_dir, "TechnologyFromStorage.csv"), index=None)
        
    res_cap.to_csv(os.path.join(output_data_dir, 'ResidualCapacity.csv'), 
                           index = None)   
  
    res_cap_storage.to_csv(os.path.join(output_data_dir, 'ResidualStorageCapacity.csv'), 
                           index = None)   

    storage_level_start.to_csv(os.path.join(output_data_dir, 'StorageLevelStart.csv'), 
                           index = None)      
    
    if tech_capacity_storage is not None:
        min_cap_invest_storage.to_csv(os.path.join(output_data_dir, 
                                                'TotalAnnualMinCapacityInvestment.csv'),
                                            index = None)
        
    else:
        min_cap_invest_base.to_csv(os.path.join(output_data_dir, 
                                                'TotalAnnualMinCapacityInvestment.csv'),
                                            index = None)

if __name__ == "__main__":
    
    if "snakemake" in globals():
        file_storage_build_rates = snakemake.input.storage_build_rates 
        file_default_op_life = snakemake.input.default_op_life
        file_gesdb_project_data = snakemake.input.gesdb_project_data
        file_gesdb_regional_mapping = snakemake.input.gesdb_regional_mapping
        storage_existing = snakemake.params.storage_existing
        storage_planned = snakemake.params.storage_planned        
        start_year = snakemake.params.start_year
        end_year = snakemake.params.end_year
        region_name = snakemake.params.region_name
        custom_nodes = snakemake.params.custom_nodes
        tech_capacity_sto = snakemake.params.user_defined_capacity_storage
        no_investment_techs = snakemake.params.no_investment_techs      
        storage_parameters = snakemake.params.storage_parameters           
        output_data_dir = snakemake.params.output_data_dir
        transmission_data_dir = snakemake.params.transmission_data_dir            
        file_iar_base = f'{transmission_data_dir}/InputActivityRatio.csv'
        file_oar_base = f'{transmission_data_dir}/OutputActivityRatio.csv'
        file_capact_base = f'{transmission_data_dir}/CapacityToActivityUnit.csv' 
        file_cap_cost_base = f'{transmission_data_dir}/CapitalCost.csv'             
        file_fix_cost_base = f'{transmission_data_dir}/FixedCost.csv'
        file_var_cost_base = f'{transmission_data_dir}/VariableCost.csv'   
        file_op_life_base = f'{transmission_data_dir}/OperationalLife.csv'          
        file_max_cap_invest_base = f'{transmission_data_dir}/TotalAnnualMaxCapacityInvestment.csv'
        file_min_cap_invest_base = f'{transmission_data_dir}/TotalAnnualMinCapacityInvestment.csv'
        file_res_cap_base = f'{transmission_data_dir}/ResidualCapacity.csv'
        file_tech_set = f'{transmission_data_dir}/TECHNOLOGY.csv'        
        file_fuel_set = f'{output_data_dir}/FUEL.csv'
        
    # The below else statement defines variables if the 'transmission/main' script is to be run locally
    # outside the snakemake workflow. This is relevant for testing purposes only! User inputs when running 
    # the full workflow need to be defined in the config file. 

    else:    
        file_storage_build_rates = 'resources/data/custom/storage_build_rates.csv'
        file_default_op_life = 'resources/data/custom/operational_life.csv'
        file_gesdb_project_data = 'resources/data/default/GESDB_Project_Data.json'
        file_gesdb_regional_mapping = 'resources/data/custom/GESDB_region_mapping.csv'
        storage_existing = True
        storage_planned = True          
        start_year = 2021
        end_year = 2050
        region_name = 'GLOBAL'
        custom_nodes = []
        tech_capacity_sto = {'sto1': ['PWRSDSINDWE01', 2, 2010, 2025, 3, 1800, 40, 0, 87],
                             'sto2': ['PWRLDSINDNE01', 4, 1985, 2025, 2, 3400, 19, 0.5, 82],
                             'sto3': ['PWRLDSINDNE01', 1, 2015, 2025, 2, 3400, 19, 0.5, 82]}
        no_investment_techs = ["CSP", "WAV", "URN", "OTH", "WAS", 
                               "COG", "GEO", "BIO", "PET", "LDS"]
        
        storage_parameters = {'SDS': [1938, 44.25, 0, 85, 4],
                              'LDS': [3794, 20.2, 0.58, 80, 10]}

        output_data_dir = 'results/data'
        transmission_data_dir = 'results/data/transmission'        
        file_iar_base = f'{transmission_data_dir}/InputActivityRatio.csv'
        file_oar_base = f'{transmission_data_dir}/OutputActivityRatio.csv'
        file_capact_base = f'{transmission_data_dir}/CapacityToActivityUnit.csv'
        file_cap_cost_base = f'{transmission_data_dir}/CapitalCost.csv'        
        file_fix_cost_base = f'{transmission_data_dir}/FixedCost.csv'
        file_var_cost_base = f'{transmission_data_dir}/VariableCost.csv'
        file_op_life_base = f'{transmission_data_dir}/OperationalLife.csv'        
        file_max_cap_invest_base = f'{transmission_data_dir}/TotalAnnualMaxCapacityInvestment.csv'
        file_min_cap_invest_base = f'{transmission_data_dir}/TotalAnnualMinCapacityInvestment.csv'
        file_res_cap_base = f'{transmission_data_dir}/ResidualCapacity.csv'
        file_tech_set = f'{transmission_data_dir}/TECHNOLOGY.csv'
        file_fuel_set = f'{output_data_dir}/FUEL.csv'

    # SET INPUT DATA
    
    if storage_parameters:
        sto_techs = list(storage_parameters.keys())
    else:
        sto_techs = None
    
    build_rates = import_storage_build_rates(file_storage_build_rates)
    
    op_life = import_op_life(file_default_op_life)
    op_life_dict = dict(zip(list(op_life['tech']),
                            list(op_life['years'])))
    
    gesdb_project_data = import_GESDB_project_data(file_gesdb_project_data)
    gesdb_regional_mapping = import_GESDB_regional_mapping(file_gesdb_regional_mapping)
    
    iar_base = import_iar_base(file_iar_base)
    oar_base = import_oar_base(file_oar_base)
    capact_base = import_capact_base(file_capact_base)
    cap_cost_base = import_cap_cost_base(file_cap_cost_base)    
    fix_cost_base = import_fix_cost_base(file_fix_cost_base)
    var_cost_base = import_var_cost_base(file_var_cost_base)
    op_life_base = import_op_life_base(file_op_life_base)       
    max_cap_invest_base = import_max_cap_invest_base(file_max_cap_invest_base)
    min_cap_invest_base = import_min_cap_invest_base(file_min_cap_invest_base)
    res_cap_base = import_res_cap_base(file_res_cap_base)  
    tech_set_base = import_set_base(file_tech_set)  
    fuel_set_base = import_set_base(file_fuel_set)  
    
    input_data = {
        "unique_sto_techs" : sto_techs,
        "default_op_life": op_life_dict,
        "gesdb_data": gesdb_project_data,
        "gesdb_mapping" : gesdb_regional_mapping,
        "build_rates" : build_rates,
        "iar_base" : iar_base,
        "oar_base" : oar_base,
        "capact_base" : capact_base,
        "cap_cost_base" : cap_cost_base,        
        "fix_cost_base" : fix_cost_base,
        "var_cost_base" : var_cost_base,
        "op_life_base" : op_life_base,
        "max_cap_invest_base" : max_cap_invest_base,
        "min_cap_invest_base" : min_cap_invest_base,
        "res_cap_base" : res_cap_base, 
        "tech_set_base" : tech_set_base,
        "fuel_set_base" : fuel_set_base,
        "tech_capacity_storage": tech_capacity_sto
    }
    
    # CALL MAIN
    main(**input_data)