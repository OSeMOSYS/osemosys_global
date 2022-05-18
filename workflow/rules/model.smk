import os

# REQUIRED 

configfile: 'config/config.yaml'

# OUTPUT FILES 

osemosys_files = os.listdir('resources/simplicity/data')

# RULES

rule geographic_filter:
    message:
        'Applying geographic filter...'
    input: 
        csv_files = expand('results/data/{osemosys_file}', osemosys_file = osemosys_files),
    params:
        config['scenario'],
        config['geographic_scope']
    output:
        csv_files = expand('results/{{scenario}}/data/{osemosys_file}', osemosys_file = osemosys_files),
        datapackage = 'results/{scenario}/datapackage.json'
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/{scenario}/logs/geographicFilter.log'
    # wildcard_constraints:
    #     scenario="[^a-zA-Z0-9_]"
    shell:
        'python workflow/scripts/osemosys_global/OPG_geographic_filter.py 2> {log}'

rule otoole_convert:
    message:
        'Creating data file...'
    input:
        datapackage = 'results/{scenario}/datapackage.json',
        csv_files = expand('results/{{scenario}}/data/{osemosys_file}', osemosys_file = osemosys_files),
    output:
        data_file = 'results/{scenario}/{scenario}.txt'
    conda:
        '../envs/otoole.yaml'
    log:
        log = 'results/{scenario}/logs/otoole_convert.log'
    shell:
        'otoole convert datapackage datafile {input.datapackage} {output} 2> {log}'

rule preprocess_data_file:
    message:
        'Preprocessing data file...'
    input:
        data_file = 'results/{scenario}/{scenario}.txt'
    output:
        data_file = 'results/{scenario}/PreProcessed_{scenario}.txt'
    conda:
        '../envs/data_processing.yaml'
    log:
        log = 'results/{scenario}/logs/preprocess_data_file.log'
    shell:
        'python resources/OSeMOSYS_GNU_MathProg/scripts/preprocess_data.py otoole {input} {output} 2> {log}'

rule create_lp_file:
    message:
        'Creating lp file...'
    input:
        model_file = 'resources/osemosys_fast_preprocessed.txt',
        data_file = 'results/{scenario}/PreProcessed_{scenario}.txt'
    output:
        lp_file = 'results/{scenario}/{scenario}.lp'
    log:
        log = 'results/{scenario}/logs/create_lp_file.log'
    shell:
        'glpsol -m {input.model_file} -d {input.data_file} --wlp {output.lp_file} --check 2> {log}'

rule solve_lp:
    message:
        'Solving model via {config[solver]}...'
    input:
        lp_file = 'results/{scenario}/{scenario}.lp'
    output:
        solution = 'results/{scenario}/{scenario}.sol'
    params:
        json = 'results/{scenario}/{scenario}.json',
        ilp = 'results/{scenario}/{scenario}.ilp'
    log:
        log = 'results/{scenario}/logs/solve_lp.log'
    shell: 
        '''
        if [ {config[solver]} = gurobi ]
        then
          gurobi_cl Method=2 ResultFile={output.solution} ResultFile={params.json} ResultFile={params.ilp} {input.lp_file}
        elif [ {config[solver]} = cplex ]
        then
          cplex -c "read {input.lp_file}" "optimize" "write {output.solution}"
        else
          cbc {input.lp_file} solve -sec 1500 -solu {output.solution}
        fi
        '''

rule transform_cplex:
    message:
        'Transforming cplex results...'
    input:
        solution = 'results/{scenario}/{scenario}.sol'
    output: 
        transform = 'results/{scenario}/{scenario}_transform.sol'
    log:
        log = 'results/{scenario}/logs/transform_cplex.log'
    shell: 
        'python resources/cplex_transform.py {input.solution} {output.transform} 2> {log}'

rule sort_cplex:
    message:
        'Sorting cplex results...'
    input:
        transform = 'results/{scenario}/{scenario}_transform.sol'
    output: 
        sort = 'results/{scenario}/{scenario}_sort.sol'
    log:
        log = 'results/{scenario}/logs/sort_cplex.log'
    shell: 
        'sort {input.transform} > {output.sort} 2> {log}'



# The script on the OSeMOSYS GitHub to convert CPLEX results into CBC results wasnt working
# for me. It is supposed to be splitting the line to find what each varaible is, but 
# the parsing doesnt seem to be working correctly...  

# I replaced the transform_cplex rule below with the uncommented one which uses the script 
# found at https://osemosys.readthedocs.io/en/latest/manual/Advanced%20functionalities.html

"""
rule transform_cplex:
    message:
        'Transforming cplex results...'
    input:
        solution = 'results/{scenario}/{scenario}.sol'
    output: 
        transformed = 'results/{scenario}/{scenario}_transformed.sol'
    params: 
        python_file = 'resources/OSeMOSYS_GNU_MathProg/scripts/convert_cplex_to_cbc.py',
        start_year = config['startYear'],
        end_year = config['endYear'],
    log:
        log = 'results/{scenario}/logs/transform_cplex.log'
    shell: 
        '''
        python {params.python_file} --start_year {params.start_year} --end_year {params.end_year} \
        {input.solution} {output.transformed} 2> {log}
        '''
"""