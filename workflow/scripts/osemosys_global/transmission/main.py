import pandas as pd
import os

from read import(
    import_plexos_2015,
    import_op_life,
    import_line_data,
    import_iar_base,
    import_oar_base,
    import_capact_base,
    import_cap_cost_base,
    import_fix_cost_base,
    import_op_life_base,
    import_max_cap_invest_base,
    import_min_cap_invest_base,
    import_res_cap_base,
    import_set_base
)

from constants import(
    DF_IAR_CUSTOM_VAL,
    DF_OAR_CUSTOM_VAL
    )

from activity import(
    activity_transmission,
    activity_transmission_limit,
    create_trn_dist_capacity_activity
    )

from costs import get_transmission_costs

from operational_life import set_op_life_transmission

from investment_constraints import cap_investment_constraints_trn

from user_defined_capacity import set_user_defined_capacity_trn

from sets import create_set_from_iterators, get_unique_fuels, get_unique_technologies

def main(
    plexos_prop: pd.DataFrame,
    default_op_life: pd.DataFrame,
    line_data: pd.DataFrame,
    interface_data: pd.DataFrame,
    iar_base: pd.DataFrame,
    oar_base: pd.DataFrame,
    capact_base: pd.DataFrame,
    cap_cost_base: pd.DataFrame,
    fix_cost_base: pd.DataFrame,
    op_life_base: pd.DataFrame,
    max_cap_invest_base: pd.DataFrame,
    min_cap_invest_base: pd.DataFrame,
    res_cap_base: pd.DataFrame,
    tech_set_base: pd.DataFrame,
    fuel_set_base: pd.DataFrame
):
    
    # CALL FUNCTIONS
    
    # Set activity ratios for transmission.
    iar_trn, oar_trn = activity_transmission(iar_base, oar_base, 
                                                               plexos_prop, interface_data,
                                                               start_year, end_year,
                                                               region_name)


    # Adjust activity limits if cross border trade is not allowed following user config.
    df_crossborder_final = activity_transmission_limit(cross_border_trade, oar_trn)
    
    # Set capital and fixed transmission costs.
    df_cap_cost_trn_final, df_fix_cost_trn_final = get_transmission_costs(line_data, oar_trn,
                                                                    cap_cost_base, fix_cost_base,
                                                                    start_year, end_year,
                                                                    region_name)

    # Set operational life for transmission.
    df_op_life_trn_final = set_op_life_transmission(iar_trn, oar_trn, 
                                                    op_life_dict, op_life_base, region_name)
    
    # Set annual capacity investment constraints.
    df_max_cap_invest_trn_final = cap_investment_constraints_trn(iar_trn, 
                                                                 max_cap_invest_base, 
                                                                 no_investment_techs, 
                                                                 start_year, 
                                                                 end_year, 
                                                                 region_name)

    # Alter output csv's based on user defined capacities following user config.
    if not tech_capacity_trn is None:
        (df_max_cap_invest_trn_final, 
         df_min_cap_invest_trn_final, 
         df_res_cap_trn_final, 
         iar_trn,
         oar_trn,
         df_op_life_trn_final, 
         df_cap_cost_trn_final
         ) = set_user_defined_capacity_trn(
            tech_capacity_trn, 
            op_life_dict, 
            min_cap_invest_base, 
            df_max_cap_invest_trn_final, 
            res_cap_base,
            iar_trn,
            oar_trn,
            df_op_life_trn_final,
            df_cap_cost_trn_final,
            start_year,
            end_year,
            region_name,
            DF_IAR_CUSTOM_VAL, 
            DF_OAR_CUSTOM_VAL
            )

    # get new additions to fuel and technology sets
    exising_techs = tech_set_base.VALUE.to_list()
    iar_techs = get_unique_technologies(iar_trn)
    oar_techs = get_unique_technologies(oar_trn)
    tech_set = create_set_from_iterators(exising_techs, iar_techs, oar_techs)
    
    exising_fuels = fuel_set_base.VALUE.to_list()
    iar_fuels = get_unique_fuels(iar_trn)
    oar_fuels = get_unique_fuels(oar_trn)
    fuel_set = create_set_from_iterators(exising_fuels, iar_fuels, oar_fuels)
    
    # assign capacity to activity unit to transmission + distribution techs
    cap_activity_trn = create_trn_dist_capacity_activity(iar_trn, oar_trn)
    cap_activity = pd.concat([capact_base, cap_activity_trn]).drop_duplicates(subset=["REGION", "TECHNOLOGY"], keep="last")
    
    # OUTPUT CSV's
    
    oar_trn.to_csv(os.path.join(output_data_dir, "OutputActivityRatio.csv"), index=None)
    
    iar_trn.to_csv(os.path.join(output_data_dir, "InputActivityRatio.csv"), index=None)
    
    df_crossborder_final.to_csv(os.path.join(output_data_dir,
                                            "TotalTechnologyModelPeriodActivityUpperLimit.csv"),
                                index = None)
    
    cap_activity.to_csv(os.path.join(output_data_dir, "CapacityToActivityUnit.csv"), index = None)
    
    df_cap_cost_trn_final.to_csv(os.path.join(output_data_dir, "CapitalCost.csv"), index = None)
    
    df_fix_cost_trn_final.to_csv(os.path.join(output_data_dir, "FixedCost.csv"), index = None)
    
    df_op_life_trn_final.to_csv(os.path.join(output_data_dir, "OperationalLife.csv"), index = None)
    
    df_max_cap_invest_trn_final.to_csv(os.path.join(output_data_dir, 
                                            'TotalAnnualMaxCapacityInvestment.csv'),
                                        index = None)

    tech_set.to_csv(os.path.join(output_data_dir, "TECHNOLOGY.csv"), index = None)
    fuel_set.to_csv(os.path.join(output_data_dir, "FUEL.csv"), index = None)
    
    if not tech_capacity_trn is None:
        
        df_res_cap_trn_final.to_csv(os.path.join(output_data_dir, 
                                                'ResidualCapacity.csv'),
                                            index = None)       
        
        df_min_cap_invest_trn_final.to_csv(os.path.join(output_data_dir, 
                                                'TotalAnnualMinCapacityInvestment.csv'),
                                            index = None)

if __name__ == "__main__":
    
    if "snakemake" in globals():
        file_plexos = snakemake.input.plexos
        file_default_op_life = snakemake.input.default_op_life
        file_line_data = snakemake.input.line_data
        start_year = snakemake.params.start_year
        end_year = snakemake.params.end_year
        region_name = snakemake.params.region_name
        custom_nodes = snakemake.params.custom_nodes
        tech_capacity_trn = snakemake.params.user_defined_capacity_transmission
        no_investment_techs = snakemake.params.no_investment_techs       
        cross_border_trade = snakemake.params.trade
        output_data_dir = snakemake.params.output_data_dir
        input_data_dir = snakemake.params.input_data_dir
        powerplant_data_dir = snakemake.params.powerplant_data_dir        
        file_iar_base = f'{powerplant_data_dir}/InputActivityRatio.csv'
        file_oar_base = f'{powerplant_data_dir}/OutputActivityRatio.csv'
        file_capact_base = f'{powerplant_data_dir}/CapacityToActivityUnit.csv'      
        file_cap_cost_base = f'{powerplant_data_dir}/CapitalCost.csv'
        file_fix_cost_base = f'{powerplant_data_dir}/FixedCost.csv'
        file_op_life_base = f'{powerplant_data_dir}/OperationalLife.csv'
        file_max_cap_invest_base = f'{powerplant_data_dir}/TotalAnnualMaxCapacityInvestment.csv'
        file_min_cap_invest_base = f'{powerplant_data_dir}/TotalAnnualMinCapacityInvestment.csv'
        file_res_cap_base = f'{powerplant_data_dir}/ResidualCapacity.csv'
        file_tech_set = f'{powerplant_data_dir}/TECHNOLOGY.csv'        
        file_fuel_set = f'{powerplant_data_dir}/FUEL.csv'        

    else:
        file_plexos = 'resources/data/PLEXOS_World_2015_Gold_V1.1.xlsx'
        file_default_op_life = 'resources/data/operational_life.csv'
        file_line_data = 'resources/data/Costs Line expansion.xlsx'
        start_year = 2021
        end_year = 2050
        region_name = 'GLOBAL'
        custom_nodes = ["INDWE", "INDEA", "INDNE", "INDNO", "INDSO"]
        tech_capacity_trn = {'TRNINDEAINDNE': [5, 1975, "open", 2030, 10, 861]}
        no_investment_techs = ["CSP", "WAV", "URN", "OTH", "WAS", 
                               "COG", "GEO", "BIO", "PET"]        
        cross_border_trade = True
        output_data_dir = 'results/data'
        input_data_dir = 'resources/data'
        powerplant_data_dir = 'results/data/powerplant'
        file_iar_base = f'{powerplant_data_dir}/InputActivityRatio.csv'
        file_oar_base = f'{powerplant_data_dir}/OutputActivityRatio.csv'
        file_capact_base = f'{powerplant_data_dir}/CapacityToActivityUnit.csv'
        file_cap_cost_base = f'{powerplant_data_dir}/CapitalCost.csv'
        file_fix_cost_base = f'{powerplant_data_dir}/FixedCost.csv'
        file_op_life_base = f'{powerplant_data_dir}/OperationalLife.csv'
        file_max_cap_invest_base = f'{powerplant_data_dir}/TotalAnnualMaxCapacityInvestment.csv'
        file_min_cap_invest_base = f'{powerplant_data_dir}/TotalAnnualMinCapacityInvestment.csv'
        file_res_cap_base = f'{powerplant_data_dir}/ResidualCapacity.csv'
        file_tech_set = f'{powerplant_data_dir}/TECHNOLOGY.csv'
        file_fuel_set = f'{powerplant_data_dir}/FUEL.csv'

    # SET INPUT DATA
    plexos_prop = import_plexos_2015(file_plexos, "prop")
    op_life = import_op_life(file_default_op_life)
    op_life_dict = dict(zip(list(op_life['tech']),
                            list(op_life['years'])))
        
    trn_line = import_line_data(file_line_data, "Lines")
    trn_interface = import_line_data(file_line_data, "Interface")
    iar_base = import_iar_base(file_iar_base)
    oar_base = import_oar_base(file_oar_base)
    capact_base = import_capact_base(file_capact_base)
    cap_cost_base = import_cap_cost_base(file_cap_cost_base)
    fix_cost_base = import_fix_cost_base(file_fix_cost_base)
    op_life_base = import_op_life_base(file_op_life_base)
    max_cap_invest_base = import_max_cap_invest_base(file_max_cap_invest_base)
    min_cap_invest_base = import_min_cap_invest_base(file_min_cap_invest_base)
    res_cap_base = import_res_cap_base(file_res_cap_base)  
    tech_set_base = import_set_base(file_tech_set)  
    fuel_set_base = import_set_base(file_tech_set)  
    
    input_data = {
        "plexos_prop": plexos_prop,
        "default_op_life": op_life,
        "line_data": trn_line,
        "interface_data": trn_interface,
        "iar_base" : iar_base,
        "oar_base" : oar_base,
        "capact_base" : capact_base,
        "cap_cost_base" : cap_cost_base,
        "fix_cost_base" : fix_cost_base,
        "op_life_base" : op_life_base,
        "max_cap_invest_base" : max_cap_invest_base,
        "min_cap_invest_base" : min_cap_invest_base,
        "res_cap_base" : res_cap_base, 
        "tech_set_base" : tech_set_base,
        "fuel_set_base" : fuel_set_base,  
    }
    
    # CALL MAIN
    main(**input_data)