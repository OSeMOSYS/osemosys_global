"""Creates a time series heat map from a model run. 

Arguments
---------
path_to_parameters : str
    File containing the parameters for generated sample
model_inputs : str
    File path to sample model inputs
model_outputs : str
    File path to model outputs
location_to_save : str
    File path to save results

Usage
-----
To run the script on the command line, type::

    python analyze_results.py path/to/parameters.csv path/to/inputs.txt 
        path/to/model/results.csv path/to/save/SA/results.csv

The ``parameters.csv`` CSV file should be formatted as follows::

    name,group,indexes,min_value,max_value,dist,interpolation_index,action
    CapitalCost,pvcapex,"GLOBAL,GCPSOUT0N",500,1900,unif,YEAR,interpolate
    DiscountRate,discountrate,"GLOBAL,GCIELEX0N",0.05,0.20,unif,None,fixed

The ``inputs.txt`` should be the output from SALib.sample.morris.sample

The ``model/results.csv`` must have an 'OBJECTIVE' column holding results OR
be a formatted output of an OSeMOSYS parameter 

"""

from pathlib import Path
from SALib.analyze import morris as analyze_morris
from SALib.plotting import morris as plot_morris
import numpy as np
import pandas as pd
import csv
import sys
import utils
import seaborn as sns
import matplotlib.pyplot as plt


from logging import getLogger

logger = getLogger(__name__)

def sort_results(df : pd.DataFrame, year : int) -> np.array:
    """Organizes a model variable results file for a morris analysis
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe grouped on model number and year
    year : int
        Year to sort results for 

    Returns
    -------
    Y : np.array
        Results for morris analysis
    """
    df = df.xs(year, level=('YEAR')).reset_index()
    # df['NUM'] = df['MODELRUN'].map(lambda x: int(x.split('_')[1]))
    df['NUM'] = df['MODELRUN']
    df = df.sort_values(by='NUM').reset_index(drop=True).drop(('NUM'), axis=1).set_index('MODELRUN')
    Y = df.to_numpy()
    return Y

def main(parameters: dict, X: np.array, model_results: pd.DataFrame, save_file: str):
    """Performs SA and plots results. 

    Parameters
    ----------
    parameters : Dict
        Parameters for generated sample
    X : np.array
        Input Sample
    model_results : pd.DataFrame
        Model results for a OSeMOSYS variable 
    save_file : str
        File path to save results
    """

    problem = utils.create_salib_problem(parameters)
    model_results = model_results.groupby(['MODELRUN','YEAR']).sum()

    years = model_results.index.unique(level='YEAR')
    SA_result_data = []

    for year in model_results.index.unique(level='YEAR'):
        Y = sort_results(model_results, year)
        Si = analyze_morris.analyze(problem, X, Y, print_to_console=False)
        SA_result_data.append(Si['mu_star'])
    
    SA_result_data = np.ma.concatenate([SA_result_data])
    columns = set([x['group'] for x in parameters])
    SA_results = pd.DataFrame(np.ma.getdata(SA_result_data), columns=columns, index=years).T

    # Save figure results
    title = Path(save_file).stem.capitalize()
    height = len(columns) / 2 + 1.5
    width = len(years) / 5
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(SA_results, cmap="coolwarm", ax=ax).set_title(title)
    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 8)
    fig.savefig(f'{save_file}.png', bbox_inches='tight')

if __name__ == "__main__":

    parameters_file = sys.argv[1]
    sample = sys.argv[2]
    result_file = sys.argv[3]
    save_file = str(Path(sys.argv[4]).with_suffix(''))
    with open(parameters_file, 'r') as csv_file:
        parameters = list(csv.DictReader(csv_file))

    X = np.loadtxt(sample, delimiter=',')
    results = pd.read_csv(result_file)

    main(parameters, X, results, save_file)
    
