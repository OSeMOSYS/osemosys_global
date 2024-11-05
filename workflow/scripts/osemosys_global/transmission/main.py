import pandas as pd
import os

#os.chdir(r'C:\Users\maart\Github\osemosys_global\workflow\scripts\osemosys_global\transmission')

from read import(
    import_gtd_existing,
    import_gtd_planned,
    import_gtd_mapping,
    import_centerpoints,
    import_transmission_build_rates,
    import_op_life,
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
    CUSTOM_TRN_BA_DICT_FROM, 
    CUSTOM_TRN_BA_DICT_TO,
    CUSTOM_TRN_BA_MISSING,
    SUBSEA_LINES,
    RETIREMENT_YEAR_TRANSMISSION,
    PLANNED_BUILD_YEAR_TRANSMISSION
    )

from data import(
    format_gtd_existing,
    format_gtd_planned,
    correct_gtd_data
    )

from activity import(
    set_transmission_losses,
    activity_transmission,
    activity_transmission_limit,
    create_trn_dist_capacity_activity
    )

from costs import get_transmission_costs

from operational_life import set_op_life_transmission

from investment_constraints import cap_investment_constraints_trn

from user_defined_capacity import set_user_defined_capacity_trn

from residual_capacity import res_capacity_transmission

from sets import(create_set_from_iterators, 
                 get_unique_fuels, 
                 get_unique_technologies)

#os.chdir(r'C:\Users\maart\Github\osemosys_global')

def main(
    default_op_life: dict[str, int],
    gtd_exist: pd.DataFrame,
    gtd_planned: pd.DataFrame,
    gtd_mapping: dict[str, str],
    centerpoints_mapping: list,
    build_rates: pd.DataFrame, 
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
    fuel_set_base: pd.DataFrame
):
    
    # CALL FUNCTIONS
    
    # Correct the input data from the Global Transmission Database (GTD).
    gtd_exist_corrected, gtd_planned_corrected = correct_gtd_data(gtd_exist, 
                                                                  gtd_planned, 
                                                                  gtd_mapping, 
                                                                  CUSTOM_TRN_BA_DICT_FROM, 
                                                                  CUSTOM_TRN_BA_DICT_TO,
                                                                  CUSTOM_TRN_BA_MISSING)
    
    # Set capital, fixed and variable transmission costs.
    cap_cost_trn, fix_cost_trn, var_cost_trn = get_transmission_costs(gtd_exist_corrected, 
                                                                      gtd_planned_corrected,
                                                                      cap_cost_base,
                                                                      fix_cost_base,
                                                                      var_cost_base,
                                                                      centerpoints_mapping, 
                                                                      transmission_parameters, 
                                                                      start_year, end_year, 
                                                                      region_name, SUBSEA_LINES)
    
    # Calculate technology and transmission pathway specific transmission losses.
    eff_trn = set_transmission_losses(gtd_exist_corrected, gtd_planned_corrected,
                                      centerpoints_mapping, transmission_parameters, SUBSEA_LINES)
    
    # Set activity ratios for transmission.
    iar_trn, oar_trn = activity_transmission(iar_base, oar_base, 
                                             eff_trn, start_year, 
                                             end_year, region_name)
    
    # Adjust activity limits if cross border trade is not allowed following user config.
    activity_limit_trn = activity_transmission_limit(cross_border_trade, oar_trn)
    
    # Set operational life for transmission.
    op_life_trn = set_op_life_transmission(oar_trn, default_op_life, op_life_base, region_name)
    
    # Set annual capacity investment constraints.
    max_cap_invest_trn = cap_investment_constraints_trn(iar_trn, 
                                                        max_cap_invest_base, 
                                                        build_rates,
                                                        no_investment_techs, 
                                                        start_year, 
                                                        end_year, 
                                                        region_name)
    
    # Set residual capacity.
    res_cap_trn = res_capacity_transmission(gtd_exist_corrected, gtd_planned_corrected, 
                                            res_cap_base, op_life_dict, 
                                            start_year, end_year, region_name,
                                            RETIREMENT_YEAR_TRANSMISSION, 
                                            PLANNED_BUILD_YEAR_TRANSMISSION)

    # Alter output csv's based on user defined capacities following user config.
    if tech_capacity_trn is not None:
        (max_cap_invest_trn, 
         min_cap_invest_trn, 
         res_cap_trn, 
         iar_trn,
         oar_trn,
         op_life_trn, 
         cap_cost_trn,
         fix_cost_trn,
         var_cost_trn
         ) = set_user_defined_capacity_trn(
            tech_capacity_trn, 
            default_op_life, 
            min_cap_invest_base, 
            max_cap_invest_trn, 
            res_cap_trn,
            iar_trn,
            oar_trn,
            op_life_trn,
            cap_cost_trn,
            fix_cost_trn,
            var_cost_trn,            
            start_year,
            end_year,
            region_name
            )  

    # get new additions to fuel and technology sets
    existing_techs = tech_set_base.VALUE.to_list()
    iar_techs = get_unique_technologies(iar_trn)
    oar_techs = get_unique_technologies(oar_trn)
    tech_set = create_set_from_iterators(existing_techs, iar_techs, oar_techs)
    
    existing_fuels = fuel_set_base.VALUE.to_list()
    iar_fuels = get_unique_fuels(iar_trn)
    oar_fuels = get_unique_fuels(oar_trn)
    fuel_set = create_set_from_iterators(existing_fuels, iar_fuels, oar_fuels)
    
    # assign capacity to activity unit to transmission + distribution techs
    cap_activity_trn = create_trn_dist_capacity_activity(iar_trn, oar_trn)
    cap_activity = pd.concat([capact_base, cap_activity_trn]
                             ).drop_duplicates(subset=["REGION", "TECHNOLOGY"], keep="last")
    
    # OUTPUT CSV's
    
    oar_trn.to_csv(os.path.join(transmission_data_dir, "OutputActivityRatio.csv"), index=None)
    
    iar_trn.to_csv(os.path.join(transmission_data_dir, "InputActivityRatio.csv"), index=None)
    
    activity_limit_trn.to_csv(os.path.join(output_data_dir,
                                            "TotalTechnologyModelPeriodActivityUpperLimit.csv"),
                                index = None)
    
    cap_activity.to_csv(os.path.join(transmission_data_dir, "CapacityToActivityUnit.csv"), index = None)
    
    cap_cost_trn.to_csv(os.path.join(output_data_dir, "CapitalCost.csv"), index = None)
    
    fix_cost_trn.to_csv(os.path.join(transmission_data_dir, "FixedCost.csv"), index = None)
    
    var_cost_trn.to_csv(os.path.join(transmission_data_dir, "VariableCost.csv"), index = None)    
    
    op_life_trn.to_csv(os.path.join(output_data_dir, "OperationalLife.csv"), index = None)
    
    max_cap_invest_trn.to_csv(os.path.join(transmission_data_dir, 
                                            'TotalAnnualMaxCapacityInvestment.csv'),
                                        index = None)

    tech_set.to_csv(os.path.join(transmission_data_dir, "TECHNOLOGY.csv"), index = None)
    fuel_set.to_csv(os.path.join(output_data_dir, "FUEL.csv"), index = None)
        
    res_cap_trn.to_csv(os.path.join(transmission_data_dir, 'ResidualCapacity.csv'),
                       index = None)       
    
    if tech_capacity_trn is not None:
        min_cap_invest_trn.to_csv(os.path.join(transmission_data_dir, 
                                                'TotalAnnualMinCapacityInvestment.csv'),
                                            index = None)
        
    else:
        min_cap_invest_base.to_csv(os.path.join(transmission_data_dir, 
                                                'TotalAnnualMinCapacityInvestment.csv'),
                                            index = None)

if __name__ == "__main__":
    
    if "snakemake" in globals():
        file_gtd_existing = snakemake.input.gtd_existing
        file_gtd_planned = snakemake.input.gtd_planned
        file_gtd_mapping = snakemake.input.gtd_mapping
        file_centerpoints = snakemake.input.centerpoints 
        file_transmission_build_rates = snakemake.input.transmission_build_rates        
        file_default_op_life = snakemake.input.default_op_life
        start_year = snakemake.params.start_year
        end_year = snakemake.params.end_year
        region_name = snakemake.params.region_name
        tech_capacity_trn = snakemake.params.user_defined_capacity_transmission
        no_investment_techs = snakemake.params.no_investment_techs      
        transmission_parameters = snakemake.params.transmission_parameters           
        cross_border_trade = snakemake.params.trade
        output_data_dir = snakemake.params.output_data_dir
        input_data_dir = snakemake.params.input_data_dir
        powerplant_data_dir = snakemake.params.powerplant_data_dir  
        transmission_data_dir = snakemake.params.transmission_data_dir 
        file_iar_base = f'{powerplant_data_dir}/InputActivityRatio.csv'
        file_oar_base = f'{powerplant_data_dir}/OutputActivityRatio.csv'
        file_capact_base = f'{powerplant_data_dir}/CapacityToActivityUnit.csv'      
        file_cap_cost_base = f'{powerplant_data_dir}/CapitalCost.csv'
        file_fix_cost_base = f'{powerplant_data_dir}/FixedCost.csv'
        file_var_cost_base = f'{powerplant_data_dir}/VariableCost.csv'        
        file_op_life_base = f'{powerplant_data_dir}/OperationalLife.csv'
        file_max_cap_invest_base = f'{powerplant_data_dir}/TotalAnnualMaxCapacityInvestment.csv'
        file_min_cap_invest_base = f'{powerplant_data_dir}/TotalAnnualMinCapacityInvestment.csv'
        file_res_cap_base = f'{powerplant_data_dir}/ResidualCapacity.csv'
        file_tech_set = f'{powerplant_data_dir}/TECHNOLOGY.csv'        
        file_fuel_set = f'{powerplant_data_dir}/FUEL.csv'
        
    # The below else statement defines variables if the 'transmission/main' script is to be run locally
    # outside the snakemake workflow. This is relevant for testing purposes only! User inputs when running 
    # the full workflow need to be defined in the config file. 

    else:
        file_gtd_existing = 'resources/data/GTD_existing.csv'
        file_gtd_planned = 'resources/data/GTD_planned.csv'    
        file_gtd_mapping = 'resources/data/GTD_region_mapping.csv'  
        file_centerpoints = 'resources/data/centerpoints.csv'         
        file_transmission_build_rates = 'resources/data/transmission_build_rates.csv'         
        file_default_op_life = 'resources/data/operational_life.csv'
        start_year = 2021
        end_year = 2050
        region_name = 'GLOBAL'

        custom_nodes = []
        tech_capacity_trn = {'trn1': ['TRNINDEAINDNE', 5, 1975, 2030, 10, 350, 13, 4, 95],
                             'trn2': ['TRNINDEAINDNE', 1, 1990, 2030, 10, 350, 13, 4, 95],
                             'trn3': ['TRNINDEAINDNE', 2, 2035, 2030, 10, 350, 13, 4, 95],
                             'trn4': ['TRNINDNOINDSO', 0, 2020, 2025, 0.5, 620, 24, 4, 92]}
        
        no_investment_techs = ["CSP", "WAV", "URN", "OTH", "WAS", 
                               "COG", "GEO", "BIO", "PET"]
        transmission_parameters = {'HVAC': [779, 95400, 6.75, 0, 3.5, 4],
                                   'HVDC': [238, 297509, 3.5, 1.3, 3.5, 4],
                                   'HVDC_subsea': [295, 297509, 3.5, 1.3, 3.5, 4]}
        cross_border_trade = True
        output_data_dir = 'results/data'
        input_data_dir = 'resources/data'
        powerplant_data_dir = 'results/data/powerplant'
        transmission_data_dir = 'results/data/transmission'
        file_iar_base = f'{powerplant_data_dir}/InputActivityRatio.csv'
        file_oar_base = f'{powerplant_data_dir}/OutputActivityRatio.csv'
        file_capact_base = f'{powerplant_data_dir}/CapacityToActivityUnit.csv'
        file_cap_cost_base = f'{powerplant_data_dir}/CapitalCost.csv'
        file_fix_cost_base = f'{powerplant_data_dir}/FixedCost.csv'
        file_var_cost_base = f'{powerplant_data_dir}/VariableCost.csv'        
        file_op_life_base = f'{powerplant_data_dir}/OperationalLife.csv'
        file_max_cap_invest_base = f'{powerplant_data_dir}/TotalAnnualMaxCapacityInvestment.csv'
        file_min_cap_invest_base = f'{powerplant_data_dir}/TotalAnnualMinCapacityInvestment.csv'
        file_res_cap_base = f'{powerplant_data_dir}/ResidualCapacity.csv'
        file_tech_set = f'{powerplant_data_dir}/TECHNOLOGY.csv'
        file_fuel_set = f'{powerplant_data_dir}/FUEL.csv'

    # SET INPUT DATA
    gtd_exist = format_gtd_existing(import_gtd_existing(file_gtd_existing))
    gtd_plan = format_gtd_planned(import_gtd_planned(file_gtd_planned))
    
    gtd_mapping = import_gtd_mapping(file_gtd_mapping)  
    gtd_mapping_dict = dict(zip(gtd_mapping['gtd_region'], 
                                     gtd_mapping['region']))
    
    centerpoints = import_centerpoints(file_centerpoints)  
    centerpoints_dict = centerpoints.to_dict('records')
    
    build_rates = import_transmission_build_rates(file_transmission_build_rates)
    
    op_life = import_op_life(file_default_op_life)
    op_life_dict = dict(zip(list(op_life['tech']),
                            list(op_life['years'])))

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
        "default_op_life": op_life_dict,
        "gtd_exist" : gtd_exist,
        "gtd_planned" : gtd_plan,
        "gtd_mapping" : gtd_mapping_dict,
        "centerpoints_mapping" : centerpoints_dict,
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
    }
    
    # CALL MAIN
    main(**input_data)