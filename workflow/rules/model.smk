import os

# REQUIRED 

configfile: 'config/config.yaml'
input_dir = config['inputDir']
output_dir = config['outputDir']
scenario = config['scenario']

# OUTPUT FILES 

otoole_output_dir = directory(Path(input_dir, scenario, 'results'))

osemosys_files = os.listdir(Path(input_dir, 'simplicity/data'))

# RULES

rule GeographicFilter:
    input: 
        expand(Path(output_dir, 'data/{osemosys_file}'), osemosys_file = osemosys_files),
        'config/config.yaml'
    output:
        expand(Path(output_dir, scenario, 'data/{osemosys_file}'), osemosys_file = osemosys_files),
        Path(output_dir, scenario, 'datapackage.json')
    conda:
        'envs/geoFilter.yaml'
    log:
        'workflow/logs/geographicFilter.log'
    shell:
        'python workflow/scripts/osemosys_global/OPG_geographic_filter.py 2> {log}'

rule otooleConvert:
    input:
        datapackage = Path(output_dir, scenario, 'datapackage.json'),
        csv_files = expand(Path(output_dir, scenario, 'data/{osemosys_file}'), osemosys_file = osemosys_files),
        default_value_csv = Path(output_dir, scenario, 'data/default_values.csv')
    output:
        Path(output_dir, scenario, f'{scenario}.txt')
    conda:
        'envs/otoole.yaml'
    log:
        'workflow/logs/otoole.log'
    shell:
        'otoole convert datapackage datafile {input.datapackage} {output} 2> {log}'

rule preProcess:
    input:
        Path(output_dir, scenario, f'{scenario}.txt')
    output:
        Path(output_dir, scenario, f'PreProcessed_{scenario}.txt')
    conda:
        'envs/preProcessData.yaml'
    log:
        'workflow/logs/preprocess.log'
    shell:
        'python {input_dir}/OSeMOSYS_GNU_MathProg/scripts/preprocess_data.py otoole {input} {output} 2> {log}'

rule createLP:
    input:
        modelFile = Path(input_dir, 'osemosys_fast_preprocessed.txt'),
        dataFile = Path(output_dir, scenario, f'PreProcessed_{scenario}.txt')
    output:
        Path(output_dir, scenario, f'{scenario}.lp')
    log:
        'workflow/logs/createLP.log'
    shell:
        'glpsol -m {input.modelFile} -d {input.dataFile} --wlp {output} --check 2> {log}'

rule cbcSolve:
    input:
        Path(output_dir, scenario, f'{scenario}.lp')
    output:
        Path(output_dir, scenario, f'{scenario}.sol')
    log:
        'workflow/logs/cbcSolve.log'
    shell: 
        'cbc {input} solve -solu {output} 2> {log}'
