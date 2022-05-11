import os

# REQUIRED 

configfile: 'config/config.yaml'
input_dir = config['inputDir']
output_dir = config['outputDir']
scenario = config['scenario']

# OUTPUT FILES 

osemosys_files = os.listdir(Path(input_dir, 'simplicity/data'))

# RULES

rule geographic_filter:
    message:
        'Applying geographic filter...'
    input: 
        expand(Path(output_dir, 'data/{osemosys_file}'), osemosys_file = osemosys_files),
        'config/config.yaml'
    output:
        expand(Path(output_dir, scenario, 'data/{osemosys_file}'), osemosys_file = osemosys_files),
        Path(output_dir, scenario, 'datapackage.json')
    conda:
        '../envs/data_processing.yaml'
    log:
        'workflow/logs/geographicFilter.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_geographic_filter.py 2> {log}'

rule user_capacity:
    message: 
        'Generating user-defined capacities...'
    input:
        Path(output_dir, scenario, 'data/TECHNOLOGY.csv'),
        Path(output_dir, scenario, 'data/TotalAnnualMinCapacityInvestment.csv'),
        Path(output_dir, scenario, 'data/TotalAnnualMaxCapacityInvestment.csv'),
        'config/config.yaml'
    output:
        expand(Path(output_dir, scenario, 'data/{output_file}'), output_file = user_capacity_files),
        touch('workflow/rules/flags/user_capacity.done')
    conda:
        '../envs/data_processing.yaml'
    log:
        'workflow/logs/user_defined_capacity.log'
    shell:
        'python workflow/scripts/osemosys_global/user_defined_capacity.py 2> {log}'

rule otoole_convert:
    message:
        'Creating data file...'
    input:
        datapackage = Path(output_dir, scenario, 'datapackage.json'),
        csv_files = expand(Path(output_dir, scenario, 'data/{osemosys_file}'), osemosys_file = osemosys_files),
        default_value_csv = Path(output_dir, scenario, 'data/default_values.csv')
    output:
        Path(output_dir, scenario, f'{scenario}.txt')
    conda:
        '../envs/otoole.yaml'
    log:
        'workflow/logs/otoole_convert.log'
    shell:
        'otoole convert datapackage datafile {input.datapackage} {output} 2> {log}'

rule preprocess_data_file:
    message:
        'Preprocessing data file...'
    input:
        Path(output_dir, scenario, f'{scenario}.txt')
    output:
        Path(output_dir, scenario, f'PreProcessed_{scenario}.txt')
    conda:
        '../envs/data_processing.yaml'
    log:
        'workflow/logs/preprocess_data_file.log'
    shell:
        'python {input_dir}/OSeMOSYS_GNU_MathProg/scripts/preprocess_data.py otoole {input} {output} 2> {log}'

rule create_lp_file:
    message:
        'Creating lp file...'
    input:
        modelFile = Path(input_dir, 'osemosys_fast_preprocessed.txt'),
        dataFile = Path(output_dir, scenario, f'PreProcessed_{scenario}.txt')
    output:
        Path(output_dir, scenario, f'{scenario}.lp')
    log:
        'workflow/logs/create_lp_file.log'
    shell:
        'glpsol -m {input.modelFile} -d {input.dataFile} --wlp {output} --check 2> {log}'

rule cbc_solve:
    message:
        'Solving model...'
    input:
        Path(output_dir, scenario, f'{scenario}.lp')
    output:
        Path(output_dir, scenario, f'{scenario}.sol')
    log:
        'workflow/logs/cbc_solve.log'
    shell: 
        'cbc {input} solve -solu {output} 2> {log}'
