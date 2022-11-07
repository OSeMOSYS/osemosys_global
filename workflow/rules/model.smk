"""Rules to build the model for each scenario"""

import os
import pandas as pd
import shutil
configfile: 'config/config.yaml'

INPUT_PARAMS = pd.read_csv('resources/otoole_files.csv')['inputs'].to_list()
INPUT_CSVS = [f'{f}.csv' for f in INPUT_PARAMS]

ruleorder: copy_default_values > geographic_filter

rule geographic_filter:
    message:
        'Applying geographic filter...'
    input: 
        csvs = expand('results/data/{osemosys_file}', osemosys_file = INPUT_CSVS),
    params:
        geographic_scope = config['geographic_scope']
    output:
        csv_files = expand('results/{{scenario}}/data/{osemosys_file}', osemosys_file = INPUT_CSVS),
    conda:
        '../envs/data_processing.yaml'
    log:
        'results/{scenario}/logs/geographicFilter.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_geographic_filter.py 2> {log}'

rule copy_datapackage:
    message:
        'Copying datapackage...'
    input:
        'resources/datapackage.json'
    output:
        'results/{scenario}/datapackage.json' 
    run:
        shutil.copyfile(input[0], output[0])

rule copy_config:
    message:
        'Copying config file...'
    input:
        'resources/config.yaml'
    output:
        'results/{scenario}/config.yaml' 
    run:
        shutil.copyfile(input[0], output[0])

rule generate_datafile:
    message:
        'Creating data file...'
    input:
        csv_files = expand('results/{{scenario}}/data/{osemosys_file}', osemosys_file = INPUT_CSVS),
        datapackage = 'results/{scenario}/datapackage.json',
        config = 'results/{scenario}/config.yaml',
    output:
        data_file = 'results/{scenario}/{scenario}.txt'
    conda:
        '../envs/otoole.yaml'
    log:
        'results/{scenario}/logs/generate_datafile.log'
    shell:
        'otoole convert datapackage datafile {input.datapackage} {output} {input.config} 2> {log}'

rule preprocess_datafile:
    message:
        'Preprocessing data file...'
    input:
        data_file = 'results/{scenario}/{scenario}.txt'
    output:
        data_file = 'results/{scenario}/PreProcessed_{scenario}.txt'
    conda:
        '../envs/data_processing.yaml'
    log:
        'results/{scenario}/logs/preprocess_datafile.log'
    shell:
        'python resources/OSeMOSYS_GNU_MathProg/scripts/preprocess_data.py otoole {input} {output} 2> {log}'

rule create_lp_file:
    message:
        'Creating lp file...'
    input:
        model_file = 'resources/osemosys_fast_preprocessed.txt',
        data_file = 'results/{scenario}/PreProcessed_{scenario}.txt'
    output:
        lp_file = temp('results/{scenario}/{scenario}.lp')
    log:
        'results/{scenario}/logs/create_lp_file.log'
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
        'results/{scenario}/logs/solve_lp.log'
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
        'results/{scenario}/logs/transform_cplex.log'
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
        'results/{scenario}/logs/sort_cplex.log'
    shell: 
        'sort {input.transform} > {output.sort} 2> {log}'