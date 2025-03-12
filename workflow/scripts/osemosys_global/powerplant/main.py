import pandas as pd
import os

from read import(
    import_plexos_2015,
    import_res_limit,
    import_build_rates,
    import_weo_regions,
    import_weo_costs,
    import_op_life,
    import_naming_convention_tech,
    import_afs,
    import_custom_res_cap,
    import_custom_res_potentials,
    import_specified_annual_demand,
)

from constants import(
    DUPLICATE_TECHS,
    MODE_LIST,
    RENEWABLES_LIST,
    PW2050_TECH_DICT
    )

from data import(
    set_generator_table,
    average_efficiency,
    )

from residual_capacity import(
    res_capacity,
    add_custom_res_cap
    )

from activity import(
    activity_master_start,
    activity_output_pwr,
    activity_input_pwr,
    activity_upstream,
    activity_master_end,
    capact
    )

from costs import(
    costs_pwr,
    costs_end,
    )

from operational_life import set_op_life

from investment_constraints import(
    cap_investment_constraints,
    set_renewable_limits,
    set_build_rates,
    set_fossil_capacity_constraints
    )

from sets import(
    create_sets,
    output_sets
    )

from user_defined_capacity import set_user_defined_capacity

from availability import(
    availability_factor,
    availability_factor_custom
    )

from renewable_targets import(
    apply_re_pct_targets,
    apply_re_abs_targets,
    )

from calibration import apply_calibration

from backstop import get_backstop_data

def main(
    plexos_prop: pd.DataFrame,
    plexos_memb: pd.DataFrame,
    res_limit: pd.DataFrame,
    build_rates: pd.DataFrame,
    weo_costs: pd.DataFrame,
    weo_regions: pd.DataFrame,
    default_op_life: pd.DataFrame,
    tech_code: pd.DataFrame,
    tech_list: pd.DataFrame,      
    tech_code_dict: pd.DataFrame,
    custom_res_cap: pd.DataFrame,
    default_af_factors: pd.DataFrame,
    specified_demand: pd.DataFrame,
):
    
    # CALL FUNCTIONS
    
    # return generator_table
    gen_table = set_generator_table(plexos_prop, plexos_memb, default_op_life, 
                                    tech_code_dict, start_year, end_year)

    # Calculate average technology efficiencies.
    df_eff_node, df_eff_tech = average_efficiency(gen_table)

    # Create master table for activity ratios.
    df_ratios = activity_master_start(gen_table, DUPLICATE_TECHS, MODE_LIST,
                                      custom_nodes, start_year, end_year)
    
    # Set OutputActivitiyRatio for powerplants and set df structure for InputActivityRatio.
    df_pwr_oar_final, df_pwr_iar_base = activity_output_pwr(df_ratios, region_name)  
    
    # Set InputActivityRatio for powerplants.
    df_pwr_iar_final = activity_input_pwr(df_pwr_iar_base, RENEWABLES_LIST, df_eff_node, 
                                         df_eff_tech, region_name)
    
    # Set OutputActivitiyRatio for upstream and international technologies/fuels.   
    df_oar_upstream, df_oar_int = activity_upstream(df_pwr_iar_final, RENEWABLES_LIST)
    
    # Combine and format activity ratios to output as csv.
    df_oar_final, df_iar_final = activity_master_end(df_pwr_oar_final, df_oar_upstream, 
                                                     df_oar_int, df_pwr_iar_final, 
                                                     DUPLICATE_TECHS)
        
    # Set capital and fixed powerplant costs.
    df_costs = costs_pwr(weo_costs)
    
    # Combine and format costs data to output as csv.
    df_cap_cost_final, df_fix_cost_final = costs_end(weo_regions, df_costs, df_oar_final)
    
    # Set CapacityToActivityUnit.
    df_capact_final = capact(df_oar_final)
 
    # Set operational life for powerplant technologies.
    df_op_life = set_op_life(tech_code_dict, df_iar_final, 
                             df_oar_final, default_op_life, region_name)
        
    # Set annual capacity investment constraints.
    df_max_cap_invest, df_min_cap_invest = cap_investment_constraints(df_iar_final, 
                                                                      no_investment_techs,
                                                                      start_year,
                                                                      end_year,
                                                                      region_name)
    # Calculate residual capacity.
    df_res_cap = res_capacity(gen_table, DUPLICATE_TECHS, start_year, end_year, region_name)

    # Adds residual capacity for custom entries.
    df_res_cap, custom_techs = add_custom_res_cap(df_res_cap, custom_res_cap, 
                                                  tech_list, start_year, 
                                                  end_year, region_name)      
    
    # Creates sets for TECHNOLOGIES including custom entries.
    tech_set = create_sets('TECHNOLOGY', df_oar_final, powerplant_data_dir, custom_techs)

    # Creates sets for FUEL.
    fuel_set = create_sets('FUEL', df_oar_final, powerplant_data_dir, [])       
    
    # Creates sets for YEAR, MODE_OF_OPERATION and REGION
    years_set, mode_list_set, regions_set = output_sets(MODE_LIST, start_year,
                                                        end_year, region_name)
    
    # Alter output csv's based on user defined capacities following user config.
    if not tech_capacity is None:
        (tech_set, 
         df_max_cap_invest, 
         df_min_cap_invest, 
         df_res_cap, 
         df_iar_final,
         df_oar_final,
         fuel_set,
         df_op_life, 
         df_capact_final, 
         df_cap_cost_final
         ) = set_user_defined_capacity(
             tech_capacity, 
             default_op_life, 
             tech_set, 
             df_min_cap_invest, 
             df_max_cap_invest, 
             df_res_cap,
             df_op_life,
             df_capact_final,
             df_cap_cost_final,
             df_iar_final,
             df_oar_final,
             fuel_set,
             start_year,
             end_year,
             region_name,
             RENEWABLES_LIST
             )
    
    # Set availability factors. Occurs after set_user_defined_capacity as tech_set gets updated.
    df_af_final = availability_factor(default_af_factors, tech_set,
                                      start_year, end_year, region_name)
    
    # Set custom availability factors if defined.
    df_af_final = availability_factor_custom(df_af_final, custom_af_factors)
    
    # Set total max capacity (potential - residuals) for renewable technologies.
    df_max_capacity = set_renewable_limits(res_limit, PW2050_TECH_DICT,
                                           custom_nodes, custom_res_potentials,
                                           df_res_cap, start_year, 
                                           end_year, region_name)
    
    # Add user defined build rates to TotalAnnualMaxCapacityInvestment
    df_max_cap_invest = set_build_rates(build_rates, tech_set, df_max_cap_invest, 
                                        df_max_capacity, start_year, end_year, region_name)
    
    # Set OAR, FUEL and AccumulatedAnnualDemand based on user defined RES generation targets.
    fuel_set, df_oar_final, df_accumulated_annual_demand = apply_re_pct_targets(res_targets, 
                                                                                geographic_scope,
                                                                                remove_nodes, 
                                                                                df_oar_final,
                                                                                RENEWABLES_LIST,
                                                                                fuel_set, 
                                                                                specified_demand,
                                                                                region_name)
    
    # Set OAR, FUEL and AccumulatedAnnualDemand after calibration (minium generation).
    fuel_set, df_oar_final, df_accumulated_annual_demand = apply_calibration(calibration, 
                                                                             df_oar_final, 
                                                                             df_accumulated_annual_demand, 
                                                                             fuel_set)
    
    # Set TotalAnnualMinCapacity based on user defined RES capacity targets.
    df_min_capacity = apply_re_abs_targets(res_targets, remove_nodes, region_name)
    
    # Set user defined fossil capacity constraints.
    df_max_capacity, df_min_capacity = set_fossil_capacity_constraints(tech_set, 
                                                                       fossil_capacity_targets, 
                                                                       df_max_capacity, 
                                                                       df_min_capacity,
                                                                       df_res_cap, 
                                                                       region_name)
    
    # add backstop technologies 
    bck_techs, bck_oar, bck_capex, bck_opex, bck_capact = get_backstop_data(tech_set, years_set, region_name)
    tech_set = pd.concat([tech_set, bck_techs])
    df_oar_final = pd.concat([df_oar_final, bck_oar])
    df_cap_cost_final = pd.concat([df_cap_cost_final, bck_capex])
    df_fix_cost_final = pd.concat([df_fix_cost_final, bck_opex])
    df_capact_final = pd.concat([df_capact_final, bck_capact])
    
    # OUTPUT CSV's USED AS INPUT FOR TRANSMISSION RULE
    
    df_res_cap.to_csv(os.path.join(powerplant_data_dir, "ResidualCapacity.csv"), index=None)
    
    df_oar_final.to_csv(os.path.join(powerplant_data_dir, "OutputActivityRatio.csv"), index=None)
    
    df_iar_final.to_csv(os.path.join(powerplant_data_dir, "InputActivityRatio.csv"), index=None)

    df_cap_cost_final.to_csv(os.path.join(powerplant_data_dir, "CapitalCost.csv"), index = None)
    
    df_fix_cost_final.to_csv(os.path.join(powerplant_data_dir, "FixedCost.csv"), index = None)
    
    df_capact_final.to_csv(os.path.join(powerplant_data_dir, "CapacityToActivityUnit.csv"), index = None)
    
    df_op_life.to_csv(os.path.join(powerplant_data_dir, "OperationalLife.csv"), index = None)
    
    df_max_cap_invest.to_csv(os.path.join(powerplant_data_dir, 
                                            'TotalAnnualMaxCapacityInvestment.csv'),
                                        index = None)
    
    df_min_cap_invest.to_csv(os.path.join(powerplant_data_dir, 
                                            'TotalAnnualMinCapacityInvestment.csv'),
                                        index = None)
    
    tech_set.to_csv(os.path.join(powerplant_data_dir, "TECHNOLOGY.csv"), index = None)
    
    fuel_set.to_csv(os.path.join(powerplant_data_dir, "FUEL.csv"), index = None)
    
    # OUTPUT CSV's NOT USED AS INPUT FOR TRANMISSION RULE
    
    df_max_capacity.to_csv(os.path.join(output_data_dir, "TotalAnnualMaxCapacity.csv"), 
                           index = None)
    
    df_accumulated_annual_demand.to_csv(os.path.join(output_data_dir, 
                                                     "AccumulatedAnnualDemand.csv"), 
                                        index=None)
    
    df_min_capacity.to_csv(os.path.join(output_data_dir, "TotalAnnualMinCapacity.csv"), 
                                        index=None)
    
    df_af_final.to_csv(os.path.join(output_data_dir, 'AvailabilityFactor.csv'), index=None)
    
    years_set.to_csv(os.path.join(output_data_dir, "YEAR.csv"), index = None)
    
    mode_list_set.to_csv(os.path.join(output_data_dir, "MODE_OF_OPERATION.csv"), index = None)
    
    regions_set.to_csv(os.path.join(output_data_dir, "REGION.csv"), index = None)

if __name__ == "__main__":
    
    if "snakemake" in globals():
        file_plexos = snakemake.input.plexos
        file_res_limit = snakemake.input.res_limit
        file_build_rates = snakemake.input.build_rates
        file_default_op_life = snakemake.input.default_op_life
        file_naming_convention_tech = snakemake.input.naming_convention_tech
        file_weo_costs = snakemake.input.weo_costs
        file_weo_regions = snakemake.input.weo_regions
        file_default_af_factors = snakemake.input.default_af_factors     
        start_year = snakemake.params.start_year
        end_year = snakemake.params.end_year
        region_name = snakemake.params.region_name
        geographic_scope = snakemake.params.geographic_scope
        custom_nodes = snakemake.params.custom_nodes
        remove_nodes = snakemake.params.remove_nodes
        tech_capacity = snakemake.params.user_defined_capacity
        no_investment_techs = snakemake.params.no_investment_techs
        custom_af_factors = snakemake.params.availability_factors
        res_targets = snakemake.params.res_targets
        fossil_capacity_targets = snakemake.params.fossil_capacity_targets      
        calibration = snakemake.params.calibration
        output_data_dir = snakemake.params.output_data_dir
        input_data_dir = snakemake.params.input_data_dir
        powerplant_data_dir = snakemake.params.powerplant_data_dir  
        file_specified_annual_demand = f'{output_data_dir}/SpecifiedAnnualDemand.csv'  
        file_custom_res_cap = snakemake.input.custom_res_cap
        file_custom_res_potentials = snakemake.input.custom_res_potentials
            
    # The below else statement defines variables if the 'powerplant/main' script is to be run locally
    # outside the snakemake workflow. This is relevant for testing purposes only! User inputs when running 
    # the full workflow need to be defined in the config file. 
            
    else:
        file_plexos = 'resources/data/default/PLEXOS_World_2015_Gold_V1.1.xlsx'
        file_res_limit = 'resources/data/default/PLEXOS_World_MESSAGEix_GLOBIOM_Softlink.xlsx'
        file_build_rates = 'resources/data/default/powerplant_build_rates.csv'        
        file_default_op_life = 'resources/data/custom/operational_life.csv'
        file_naming_convention_tech = 'resources/data/default/naming_convention_tech.csv'
        file_weo_costs = 'resources/data/default/weo_2020_powerplant_costs.csv'
        file_weo_regions = 'resources/data/default/weo_region_mapping.csv'
        file_default_af_factors = 'resources/data/custom/availability_factors.csv'        
        start_year = 2021
        end_year = 2050
        region_name = 'GLOBAL'
        geographic_scope = ['BTN', 'IND']
        custom_nodes = [] 
        remove_nodes = []
        tech_capacity = {'PWRCOAINDWE01': [8, 2000, 2025, 5, 1100, 35],
                         'PWRBIOINDWE01': [0, 2020, 2030, 2, 2000, 28]}
        no_investment_techs = ["CSP", "WAV", "URN", "OTH", "WAS", 
                               "COG", "GEO", "BIO", "PET"]
        custom_af_factors =  [["INDWE", 'COA', 2023, 2050, 50], 
                              ["IND", 'COA', 2023, 2050, 25]]
        res_targets = {'T01': ["", [], "PCT", 2048, 2050, 95],
                       'T02': ["IND", [], "PCT", 2030, 2040, 60],
                       'T03': ["INDSO", ['WOF','WON'], "PCT", 2025, 2045, 15],
                       'T04': ["INDSO", ['WOF'], "ABS", 2040, 2050, 100] 
                      }
        fossil_capacity_targets = [["BTNXX", 'COA', 2030, 2050, 'ABS', 1],
                                   ["INDNE", 'CCG', 2040, 2050, 'MIN', 10],
                                   ["INDSO", 'OCG', 2025, 2050, 'MAX', 25]]
        calibration = {'OCG': [50, "IND", 2021]}
        output_data_dir = 'results/data'
        input_data_dir = 'resources/data'
        powerplant_data_dir = 'results/data/powerplant'
        file_specified_annual_demand = f'{output_data_dir}/SpecifiedAnnualDemand.csv'
        file_custom_res_cap = 'resources/data/custom/residual_capacity.csv'
        file_custom_res_potentials = 'resources/data/custom/RE_potentials.csv' 

    # SET INPUT DATA
    plexos_prop = import_plexos_2015(file_plexos, "prop")
    plexos_memb = import_plexos_2015(file_plexos, "memb")
    
    res_limit = import_res_limit(file_res_limit)
    build_rates = import_build_rates(file_build_rates)
    
    op_life = import_op_life(file_default_op_life)
    
    op_life_dict = dict(zip(list(op_life['tech']),
                            list(op_life['years'])))
        
    tech_code = import_naming_convention_tech(file_naming_convention_tech)      
    
    tech_list = list(tech_code['code'])
    
    tech_code_dict = dict(zip(list(tech_code['tech']),
                              list(tech_code['code'])))
    
    weo_costs = import_weo_costs(file_weo_costs)
    weo_regions = import_weo_regions(file_weo_regions)
    
    availability = import_afs(file_default_af_factors)
    specified_annual_demand = import_specified_annual_demand(file_specified_annual_demand)

    custom_res_cap = import_custom_res_cap(file_custom_res_cap)
    custom_res_potentials = import_custom_res_potentials(file_custom_res_potentials)
    
    input_data = {
        "plexos_prop": plexos_prop,
        "plexos_memb": plexos_memb,
        "res_limit" : res_limit,
        "build_rates" : build_rates,
        "weo_costs": weo_costs,
        "weo_regions": weo_regions,
        "default_op_life": op_life_dict,
        "tech_code": tech_code,
        "tech_list": tech_list,        
        "tech_code_dict": tech_code_dict,
        "custom_res_cap" : custom_res_cap,
        "default_af_factors": availability,
        "specified_demand" : specified_annual_demand
    }
    
    # CALL MAIN
    main(**input_data)