"""Rules for Sensitivity Analysis"""

GROUPS = pd.read_csv(config['sensitivity_analysis']['parameters'])['group'].unique()
# Calculates number of model runs for the Method of Morris
MODELRUNS = range((len(GROUPS) + 1) * config['sensitivity_analysis']['replicates'])
SA_RESULTS = pd.read_csv(config['sensitivity_analysis']["results"])
SA_RESULT_FILES = SA_RESULTS['filename'].to_list()
ZIP = '.gz' if config['sensitivity_analysis']['zip'] else ''

def get_SA_input(wildcards):
    input_file = SA_RESULTS.set_index('filename').loc[wildcards.result_file]['resultfile']
    return ["results/{wildcards.scenario}/SA/modelruns/model_{modelrun}/results/{input_file}.csv".format(
        modelrun=x, input_file=input_file, wildcards=wildcards) for x in MODELRUNS]

def get_indices(wildcards):
    indices = SA_RESULTS.set_index('filename').loc[wildcards.result_file].dropna().drop('resultfile').to_dict()
    return {x:str(indices[x]).split(',') for x in indices}

def solver_file_type_sa(wildcards):
    if config['solver'] == 'cplex':
        return 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}_sort.sol'
    else: 
        return 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.sol'

rule create_sample:
    message: "Creating sample for '{params.replicates}' trajectories and '{params.parameters}' parameters"
    params:
        replicates=config['sensitivity_analysis']['replicates'],
        parameters=config['sensitivity_analysis']['parameters']
    output: "results/{scenario}/SA/modelruns/morris_sample.txt"
    # conda: "envs/sample.yaml"
    log: "results/log/create_{scenario}_sample.log"
    shell:
        "python workflow/scripts/osemosys_global/sensitivity/create_sample.py {params.parameters} {output} {params.replicates}"

rule expand_sample:
    params:
        parameters=config['sensitivity_analysis']['parameters']
    input: "results/{scenario}/SA/modelruns/morris_sample.txt"
    output: expand("results/{{scenario}}/SA/modelruns/model_{model_run}/sample_{model_run}.txt", model_run=MODELRUNS)
    # conda: "envs/sample.yaml"
    log: "results/log/expand_{scenario}_sample.log"
    shell:
        "python workflow/scripts/osemosys_global/sensitivity/expand_sample.py {input} {params.parameters} {output}"

rule create_model_data:
    message: "Copying and modifying data for '{params.folder}'"
    input:
        csvs=expand('results/{{scenario}}/data/{osemosys_file}', osemosys_file = INPUT_CSVS),
        datapackage='results/{scenario}/datapackage.json',
        sample="results/{scenario}/SA/modelruns/model_{model_run}/sample_{model_run}.txt",
        config='results/{scenario}/config.yaml',
    log: "results/log/copy_datapackage_{scenario}_{model_run}.log"
    params:
        folder=directory("results/{scenario}/SA/modelruns/model_{model_run}")
    # conda: "../envs/otoole.yaml"
    group: "gen_lp"
    output:
        csv=expand("results/{{scenario}}/SA/modelruns/model_{{model_run}}/data/{csv}.csv", csv=INPUT_FILES)
    shell:
        "python workflow/scripts/osemosys_global/sensitivity/create_modelrun.py {input.datapackage} {params.folder}/data {input.sample} {input.config}"

rule generate_SA_datafile:
    message:
        'Creating data file...'
    input:
        csv_files = expand("results/{{scenario}}/SA/modelruns/model_{{model_run}}/data/{csv}.csv", csv = INPUT_FILES),
        datapackage = 'results/{scenario}/datapackage.json',
        config = 'results/{scenario}/config.yaml',
    output:
        data_file = "results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.txt"
    # conda:
    #     '../envs/otoole.yaml'
    log:
        'results/{scenario}/logs/{model_run}/generate_datafile.log'
    shell:
        'otoole convert datapackage datafile {input.datapackage} {output} {input.config} 2> {log}'

rule preprocess_SA_datafile:
    message:
        'Preprocessing data file...'
    input:
        data_file = "results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.txt"
    output:
        data_file = "results/{scenario}/SA/modelruns/model_{model_run}/PreProcessed_{model_run}.txt"
    # conda:
    #     '../envs/data_processing.yaml'
    log:
        'results/{scenario}/logs/preprocess_datafile_{scenario}_{model_run}.log'
    shell:
        'python resources/OSeMOSYS_GNU_MathProg/scripts/preprocess_data.py otoole {input} {output} 2> {log}'

rule create_SA_lp_file:
    message:
        'Creating lp file...'
    input:
        model_file = 'resources/osemosys_fast_preprocessed.txt',
        data_file = "results/{scenario}/SA/modelruns/model_{model_run}/PreProcessed_{model_run}.txt"
    output:
        lp_file = temp('results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.lp')
    log:
        'results/{scenario}/logs/create_lp_file_{scenario}_{model_run}.log'
    shell:
        'glpsol -m {input.model_file} -d {input.data_file} --wlp {output.lp_file} --check 2> {log}'

rule solve_SA_lp:
    message:
        'Solving model via {config[solver]}...'
    input:
        lp_file = 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.lp'
    output:
        solution = 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.sol'
    params:
        json = 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.json',
        ilp = 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.ilp'
    log:
        'results/{scenario}/logs/solve_lp_{scenario}_{model_run}.log'
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

rule transform_SA_cplex:
    message:
        'Transforming cplex results...'
    input:
        solution = 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.sol'
    output: 
        transform = 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}_transform.sol'
    log:
        'results/{scenario}/logs/transform_cplex_{scenario}_{model_run}.log'
    shell: 
        'python resources/cplex_transform.py {input.solution} {output.transform} 2> {log}'

rule sort_SA_cplex:
    message:
        'Sorting cplex results...'
    input:
        transform = 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}_transform.sol'
    output: 
        sort = 'results/{scenario}/SA/modelruns/model_{model_run}/{model_run}_sort.sol'
    log:
        'results/{scenario}/logs/sort_cplex_{scenario}_{model_run}.log'
    shell: 
        'sort {input.transform} > {output.sort} 2> {log}'

rule unzip:
    message: "Unzipping LP file"
    input:
        "results/{scenario}/SA/temp/model_{model_run}.lp.gz"
    group:
        "solve"
    output:
        temp("results/{scenario}/SA/temp/model_{model_run}.lp")
    shell:
        "gunzip -fcq {input} > {output}"

rule zip_solution:
    message: "Zip up solution file {input}"
    group: "solve"
    input: "results/{scenario}/SA/temp/model_{model_run}.sol"
    output: expand("{{scenario}}/results/SA/temp/{{model_run}}.sol{zip_extension}", zip_extension=ZIP)
    shell: "gzip -fcq {input} > {output}"

rule unzip_solution:
    message: "Unzip solution file {input}"
    group: "results"
    input: "results/{scenario}/SA/temp/{model_run}.sol.gz"
    output: temp("results/{scenario}/SA/temp/{model_run}.sol")
    shell: "gunzip -fcq {input} > {output}"

rule get_statistics:
    message: "Extract the {config[solver]} statistics from the sol file"
    input: rules.solve_SA_lp.output.solution
    output: "results/{scenario}/SA/modelruns/model_{model_run}/{model_run}.stats"
    group: "solve"
    shell: 
        """
        if [ {config[solver]} = cplex ]
        then
          head -n 27 {input} | tail -n 25 > {output}
        else
          head -n 1 {input} > {output}
        fi
        """

rule get_objective_value:
    input: expand("results/{{scenario}}/SA/modelruns/model_{model_run}/{model_run}.stats",  model_run=MODELRUNS)
    output: "results/{scenario}/SA/objective.csv"
    shell:
        """
        echo "FILE,OBJECTIVE,STATUS" > {output}
        if [ {config[solver]} = cplex ]
        then
          for FILE in {input}
          do
          OBJ=$(head $FILE | grep -e 'objectiveValue' | cut -f 2 -d '=')
          STATUS=$(head $FILE | grep -e 'solutionStatusString' | cut -f 2 -d '=')
          JOB=$(echo $FILE | cut -f 3 -d '/' | cut -f 1 -d '.')
          echo "$JOB,$OBJ,$STATUS" >> {output}
          done
        elif [ {config[solver]} = cbc ]
        then
          i=0
          for FILE in {input}
          do
          OBJ=$(head $FILE | cut -f 5 -d ' ')
          STATUS=$(head $FILE | cut -f 1 -d ' ')
          JOB=$FILE
          echo "$JOB,$OBJ,$STATUS" >> {output}
          ((i=i+1))
          done
        else
          echo "To be done"
        fi
        """

rule otoole_SA_results:
    message:
        'Generating result csv files...'
    input:
        solution_file = solver_file_type_sa,
        data_file = "results/{scenario}/SA/modelruns/model_{model_run}/PreProcessed_{model_run}.txt"
    params: 
        datapackage = 'results/{scenario}/datapackage.json',
        config = 'results/{scenario}/config.yaml',
        folder = 'results/{scenario}/SA/modelruns/model_{model_run}/results'
    output:
        expand("results/{{scenario}}/SA/modelruns/model_{{model_run}}/results/{result}.csv", result = OUTPUT_FILES),
    # conda:
    #     '../envs/otoole.yaml'
    log:
        log = 'results/{scenario}/logs/otoole_results_{scenario}_{model_run}.log'
    shell: 
        '''
        otoole results {config[solver]} csv \
        {input.solution_file} {params.folder} \
        --input_datafile {input.data_file} \
        --input_datapackage {params.datapackage} \
        {params.config} 2> {log} 
        '''

rule extract_results:
    input: 
        csvs=get_SA_input,
        config='results/{scenario}/config.yaml',
    params:
        parameter = get_indices,
        folder=directory("results/{{scenario}}/SA/results/")
    log: "results/log/extract_scenario_{scenario}_{result_file}.log"
    output: expand("results/{{scenario}}/SA/modelruns/{{result_file}}.{ext}", ext=config['sensitivity_analysis']['filetype'])
    # conda: "../envs/otoole.yaml"
    script: "../scripts/osemosys_global/sensitivity/extract_results.py"

rule calculate_SA_objective:
    message:
        "Calcualting objective cost sensitivity measures"
    params: 
        parameters=config['sensitivity_analysis']['parameters'],
        result_type='objective'
    input: 
        sample = "results/{scenario}/SA/modelruns/morris_sample.txt",
        results = "results/{scenario}/SA/objective.csv"
    output: 
        expand("results/{{scenario}}/SA/objective_results.{ext}",ext=['csv','png'])
    # conda: "../envs/sample.yaml"
    shell: "python workflow/scripts/osemosys_global/sensitivity/calculate_SA_results.py {params.parameters} {input.sample} {input.results} {output[0]} {params.result_type}"

rule calculate_SA_user_defined:
    message:
        "Calcualting user defined sensitivity measures"
    params: 
        parameters=config['sensitivity_analysis']['parameters'],
        result_type='variable'
    input: 
        sample = "results/{scenario}/SA/modelruns/morris_sample.txt",
        results=expand("results/{{scenario}}/SA/modelruns/{{result_file}}.{ext}", ext=config['sensitivity_analysis']['filetype'])
    output: 
        expand("results/{{scenario}}/SA/user_results/{{result_file}}.{ext}",ext=['csv','png'])
    # conda: "../envs/sample.yaml"
    shell: "python workflow/scripts/osemosys_global/sensitivity/calculate_SA_results.py {params.parameters} {input.sample} {input.results} {output[0]} {params.result_type}"

rule create_heatmap:
    message: 
        "Calculating user defined sensitivity measures"
    params: 
        parameters=config['sensitivity_analysis']['parameters']
    input:
        sample="results/{scenario}/SA/modelruns/morris_sample.txt",
        results=expand("results/{{scenario}}/SA/modelruns/{{result_file}}.{ext}", ext=config['sensitivity_analysis']['filetype'])
    output:
        "results/{scenario}/SA/user_results/{result_file}_heatmap.png"
    shell: "python workflow/scripts/osemosys_global/sensitivity/create_heatmap.py {params.parameters} {input.sample} {input.results} {output}"