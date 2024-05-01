import os
import shutil

# REQUIRED 

configfile: 'config/config.yaml'

# OUTPUT FILES 

#osemosys_files = os.listdir('resources/simplicity/data')
osemosys_files = os.listdir('resources/otoole/data')

# RULES

rule geographic_filter:
    message:
        'Applying geographic filter...'
    input: 
        csv_files = expand('results/data/{osemosys_file}', osemosys_file = osemosys_files),
    params:
        geographic_scope = config['geographic_scope']
    output:
        csv_files = expand('results/{{scenario}}/data/{osemosys_file}', osemosys_file = osemosys_files),
        # datapackage = 'results/{scenario}/datapackage.json'
    # conda:
    #     '../envs/data_processing.yaml'
    log:
        log = 'results/{scenario}/logs/geographicFilter.log'
    shell:
        'python workflow/scripts/osemosys_global/geographic_filter.py 2> {log}'

rule copy_otoole_confg:
    message:
        'Copying otoole configuration file...'
    input:
        config='resources/otoole/config.yaml'
    output:
        config='results/{scenario}/otoole.yaml'
    run:
        shutil.copyfile(input.config, output.config)

rule copy_og_config:
    message:
        'Copying OSeMOSYS Global Configuration File'
    input:
        config='config/config.yaml'
    output:
        config='results/{scenario}/og.yaml'
    run:
        shutil.copyfile(input.config, output.config)

rule otoole_convert:
    message:
        'Creating data file...'
    params:
        csv_dir = 'results/{scenario}/data/'
    input:
        otoole_config = 'results/{scenario}/otoole.yaml',
        csv_files = expand('results/{{scenario}}/data/{osemosys_file}', osemosys_file = osemosys_files),
    output:
        data_file = 'results/{scenario}/{scenario}.txt'
    log:
        log = 'results/{scenario}/logs/otoole_convert.log'
    shell:
        'otoole convert csv datafile {params.csv_dir} {output} {input.otoole_config} 2> {log}'

rule preprocess_data_file:
    message:
        'Preprocessing data file...'
    input:
        data_file = 'results/{scenario}/{scenario}.txt'
    output:
        data_file = 'results/{scenario}/PreProcessed_{scenario}.txt'
    #conda:
    #    '../envs/data_processing.yaml'
    log:
        log = 'results/{scenario}/logs/preprocess_data_file.log'
    shell:
        'python resources/preprocess_data.py {input} {output} 2> {log}'

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
        solution = 'results/{scenario}/{scenario}.sol',
    params:
        json = 'results/{scenario}/{scenario}.json',
        ilp = 'results/{scenario}/{scenario}.ilp',
        duals = 'results/{scenario}/{scenario}.attr'
    log:
        log = 'results/{scenario}/logs/solve_lp.log'
    shell: 
        '''
        if [ {config[solver]} = gurobi ]
        then
          gurobi_cl Method=2 ResultFile={output.solution} ResultFile={params.duals} ResultFile={params.json} ResultFile={params.ilp} {input.lp_file}
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