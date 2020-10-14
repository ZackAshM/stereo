"""
Plots the data obtained from tetra3_sweep.py
"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
import argparse
import os.path as op
data_dir = op.join(op.dirname(op.dirname(op.abspath(__file__))), *("data", "sim"))

# plot histogram of indices
def param_hist(index, filename=None):
    """index=the parameter index column, filename=save file name if not none"""
    
    len_index = int(np.max(index) - np.min(index) + 1)

    fig, ax = plt.subplots(figsize=[20, 14])
    ax.set(xlabel='Parameter Index', ylabel='Counts', 
           title=r"Best Parameter Set by Count")
    ax.grid(True)

    # the histogram of the data
    n, bins, patches = plt.hist(index, len_index, density=False, facecolor='darkgreen', alpha=0.8, rwidth=0.8)
    
    if filename:
        plt.savefig(filename+'_ParamHist.png')

    #plot
    # plt.show()
    
    # return the indices with the max solve counts
    return np.array([i for i in range(len(n)) if n[i]==np.max(n)])


# get the best parameter set via min solve time
def get_best_params(data, filename=None):
    
    len_ind = int(np.max(data[0]) - np.min(data[0]) + 1)
    
    # get the indices with max solve counts
    maxct_ind = param_hist(data[0], filename=filename)
    
    # get the params corresponding to maxct_ind
    maxct_data = data[:,np.isin(data[0],maxct_ind)]
    
    # plot time vs ind
    solve_time = np.array([np.mean(maxct_data[6][maxct_data[0]==ind]) for ind in maxct_ind])
    solve_time_err = np.array([np.std(maxct_data[6][maxct_data[0]==ind]) for ind in maxct_ind])
    
    plot_data = np.array([[x, y] for x,y in sorted(zip(maxct_ind, solve_time))], dtype=float).T
    
    fig, ax = plt.subplots(figsize=[15, 10])

    ax.set(xlabel='Parameter Index', ylabel='Solve Time [ms]', 
           title='tetra3 Solve Time vs Parameter Set')

    plt.errorbar(*plot_data, yerr=solve_time_err, fmt='o', capsize=5, c='blue')
    
    if filename:
        plt.savefig(filename+'_time.png')
    
    # plt.show()
    
    # get the min
    best_data = plot_data.T[plot_data[1] == plot_data[1].min()]
    
    best_param_ind = best_data[:,0]
    best_time = best_data[:,1]
    best_time_err = solve_time_err[plot_data[1] == plot_data[1].min()]
    
    # fov, pme, pcs, mr, mt
    best_params = np.array([[data[t][data[0]==ind][0] for ind in best_param_ind] for t in range(1, 6)])
    
    if filename:
        print(filename + ' SWEEP PARAMS:\n' + 
              'ind = {0}\nfov = {1}\npme = {2}\npcs = {3}\nmr = {4}\nmt = {5}\ntime = {6} +/- {7}'.format(
                  best_param_ind, *best_params, best_time, best_time_err)
             )
    else:
        print('SWEEP PARAMS:\n' +
              'ind = {0}\nfov = {1}\npme = {2}\npcs = {3}\nmr = {4}\nmt = {5}\ntime = {6} +/- {7}'.format(
                  best_param_ind, *best_params, best_time, best_time_err)
             )
    
    return best_params, best_time

# parse the command line arguments
parser = argparse.ArgumentParser(
    description='Plots the solve count histogram across the different parameter sets as well as the solve time across the sets with the most solved images. Prints the parameter set with the minimum solve time.')
parser.add_argument('DATA_FILE', type=str,
                    help='str; the name of the data file, .txt extension assumed')
parser.add_argument('SAVE_NAME', type=str,
                    help='str; the base name for the saved pngs')
args = parser.parse_args()

FILE = "{0}/{1}.txt".format(data_dir, args.DATA_FILE)
data = np.genfromtxt(FILE, unpack=True)

SAVE = "{0}/{1}".format(data_dir, args.SAVE_NAME)

x = get_best_params(data, filename=SAVE)

print(x)
