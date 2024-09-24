"""Creates demand projections"""

import pandas as pd
import os

from read import(
    import_plexos_2015,
    import_weo_regions,
    import_weo_costs,
    import_op_life,
    import_naming_convention_tech,
    import_line_data,
    import_afs,
    import_custom_res_cap,
)

from constants import(
    file_plexos,
    file_default_op_life,
    file_naming_convention_tech,
    file_weo_costs,
    file_weo_regions,
    file_line_data,
    file_custom_res_cap,
    file_default_av_factors,  
    avg_csp_eff,
    avg_urn_eff,
    duplicate_techs,
    nodes_extra_list,
    mode_list,
    thermal_fuel_list_oar,
    thermal_fuel_list_iar,
    renewables_list,
    thermal_fuel_list_mining,
    costs_dict,
    cross_border_trade,
    no_investment_techs,
    tech_capacity,
    df_iar_custom_val,
    df_oar_custom_val,
    custom_nodes,
    output_data_dir,
    )

from data import(
    set_generator_table,
    average_efficiency,
    )

from residual_capacity import res_capacity

from activity import(
    activity_master_start,
    activity_output_pwr,
    activity_input_pwr,
    activity_upstream,
    activity_master_end,
    capact
    )

from activity_transmission import(
    activity_transmission,
    activity_transmission_limit
    )

from costs import(
    costs_pwr,
    costs_end
    )

from costs_transmission import get_transmission_costs

from operational_life import(
    set_op_life,
    set_op_life_transmission
    )

from investment_constraints import cap_investment_constraints

from sets import(
    create_sets,
    output_sets
    )

from user_defined_capacity import set_user_defined_capacity

from availability import availability_factor

def main(
    plexos_prop: pd.DataFrame,
    plexos_memb: pd.DataFrame,
    weo_costs: pd.DataFrame,
    weo_regions: pd.DataFrame,
    default_op_life: pd.DataFrame,
    naming_convention_tech: pd.DataFrame,
    line_data: pd.DataFrame,
    interface_data: pd.DataFrame,
    custom_res_cap: pd.DataFrame,
    default_av_factors: pd.DataFrame,
):
    
    # CALL FUNCTIONS
    
    # return generator_table
    gen_table = set_generator_table(plexos_prop, plexos_memb, 
                                op_life_dict, tech_code_dict)

    # Calculate average technology efficiencies.
    df_eff_node, df_eff_tech = average_efficiency(gen_table, avg_csp_eff, 
                                                  avg_urn_eff)
    
    # Calculate residual capacity.
    df_res_cap, custom_techs = res_capacity(gen_table, tech_list, tech_code, 
                                                 custom_res_cap, duplicate_techs)

    df_ratios = activity_master_start(gen_table, nodes_extra_list, 
                                      duplicate_techs, mode_list)
    
    df_oar, df_oar_base = activity_output_pwr(df_ratios, thermal_fuel_list_oar)  
    
    df_iar_base = activity_input_pwr(df_oar, thermal_fuel_list_iar, renewables_list, 
                                     df_eff_node, df_eff_tech)
    
    df_oar_upstream, df_oar_int = activity_upstream(df_iar_base, renewables_list, 
                                                    thermal_fuel_list_mining)

    df_iar_trn, df_oar_trn, df_int_trn_oar, df_int_trn_iar = activity_transmission(df_oar_base, 
                                                                                   plexos_prop, 
                                                                                   interface_data)
    
    df_oar_final, df_iar_final = activity_master_end(df_oar_base, df_oar_upstream, df_oar_int, 
                                                     df_oar_trn, df_int_trn_oar, df_iar_base, 
                                                     df_iar_trn, df_int_trn_iar, duplicate_techs)
    
    df_costs = costs_pwr(weo_costs, costs_dict)
    
    df_trans_capex, df_trans_fix = get_transmission_costs(trn_line, df_oar_final)

    df_cap_cost_final, df_fix_cost_final = costs_end(weo_regions, df_costs, df_oar_final, 
                                             df_trans_capex, df_trans_fix)
    
    df_capact_final = capact(df_oar_final)
    
    df_crossborder_final = activity_transmission_limit(cross_border_trade, df_oar_final)
 
    df_op_life_trn, df_op_life_out = set_op_life(tech_code_dict, df_iar_final, 
                                       df_oar_final, op_life_dict)
    
    df_op_life = set_op_life_transmission(df_op_life_trn, df_op_life_out, op_life_dict)
    
    df_max_cap_invest, df_min_cap_invest = cap_investment_constraints(df_iar_final, 
                                                                      no_investment_techs)
    
    # Creates sets for TECHNOLOGIES and FUELS.
    if custom_nodes:
        tech_set = create_sets('TECHNOLOGY', df_oar_final, output_data_dir, custom_techs)
    else:
        tech_set = create_sets('TECHNOLOGY', df_oar_final, output_data_dir, [])
        
    fuel_set = create_sets('FUEL', df_oar_final, output_data_dir, [])       
    
    # Create sets for YEAR, MODE_OF_OPERATION and REGION
    years_set, mode_list_set, regions_set = output_sets(mode_list)
    
    if not tech_capacity is None:
        (tech_set, 
         fuel_set, 
         df_max_cap_invest, 
         df_min_cap_invest, 
         df_res_cap, 
         df_iar_final, 
         df_oar_final, 
         df_op_life, 
         df_capact_final, 
         df_cap_cost_final
         ) = set_user_defined_capacity(
             tech_capacity, 
             op_life_dict, 
             tech_set, 
             df_min_cap_invest, 
             df_max_cap_invest, 
             df_res_cap,
             df_iar_final,
             df_oar_final,
             fuel_set,
             df_op_life,
             df_capact_final,
             df_cap_cost_final, 
             df_iar_custom_val, 
             df_oar_custom_val
             )
    
    df_af_final = availability_factor(availability, tech_set)

    # OUTPUT CSV's
    
    df_res_cap.to_csv(os.path.join(output_data_dir, "ResidualCapacity.csv"), index=None)
    
    df_oar_final.to_csv(os.path.join(output_data_dir, "OutputActivityRatio.csv"), index=None)
    
    df_iar_final.to_csv(os.path.join(output_data_dir, "InputActivityRatio.csv"), index=None)

    df_cap_cost_final.to_csv(os.path.join(output_data_dir, "CapitalCost.csv"), index = None)
    
    df_fix_cost_final.to_csv(os.path.join(output_data_dir, "FixedCost.csv"), index = None)
    
    df_capact_final.to_csv(os.path.join(output_data_dir, "CapacityToActivityUnit.csv"), index = None)
    
    df_crossborder_final.to_csv(os.path.join(output_data_dir,
                                            "TotalTechnologyModelPeriodActivityUpperLimit.csv"),
                                index = None)
    
    df_op_life.to_csv(os.path.join(output_data_dir, "OperationalLife.csv"), index = None)
    
    df_max_cap_invest.to_csv(os.path.join(output_data_dir, 
                                            'TotalAnnualMaxCapacityInvestment.csv'),
                                        index = None)
    
    df_min_cap_invest.to_csv(os.path.join(output_data_dir, 
                                            'TotalAnnualMinCapacityInvestment.csv'),
                                        index = None)
    
    df_af_final.to_csv(os.path.join(output_data_dir, 'AvailabilityFactor.csv'), index=None)
    
    tech_set.to_csv(os.path.join(output_data_dir, "TECHNOLOGY.csv"), index = None)
    
    fuel_set.to_csv(os.path.join(output_data_dir, "FUEL.csv"), index = None)
    
    years_set.to_csv(os.path.join(output_data_dir, "YEAR.csv"), index = None)
    
    mode_list_set.to_csv(os.path.join(output_data_dir, "MODE_OF_OPERATION.csv"), index = None)
    
    regions_set.to_csv(os.path.join(output_data_dir, "REGION.csv"), index = None)

if __name__ == "__main__":

    plexos_prop = import_plexos_2015(file_plexos, "prop")
    plexos_memb = import_plexos_2015(file_plexos, "memb")
    op_life = import_op_life(file_default_op_life)
    
    op_life_dict = dict(zip(list(op_life['tech']),
                            list(op_life['years'])))
        
    tech_code = import_naming_convention_tech(file_naming_convention_tech)      
    
    tech_list = list(tech_code['code'])
    
    tech_code_dict = dict(zip(list(tech_code['tech']),
                              list(tech_code['code'])))
    
    trn_line = import_line_data(file_line_data, "Lines")
    trn_interface = import_line_data(file_line_data, "Interface")
    
    weo_costs = import_weo_costs(file_weo_costs)
    weo_regions = import_weo_regions(file_weo_regions)
    
    custom_res_cap = import_custom_res_cap(file_custom_res_cap)
    
    availability = import_afs(file_default_av_factors)
    
    input_data = {
        "plexos_prop": plexos_prop,
        "plexos_memb": plexos_memb,
        "weo_costs": weo_costs,
        "weo_regions": weo_regions,
        "default_op_life": op_life,
        "naming_convention_tech": tech_code,
        "line_data": trn_line,
        "interface_data": trn_interface,
        "custom_res_cap" : custom_res_cap,
        "default_av_factors": availability,
    }

    main(**input_data)