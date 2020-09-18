"""
Plots the data obtained from run_sim.py or run_sim.sh
"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
import argparse
import os.path as op

def analyze_alt_error(alt_err, bins='auto', arcsec=True, save=None):
    """
    alt_err=data column of altitude errors, bins=hist bins, arcsec=if True, plot alt_err 
    in arcsec instead of degrees, save=if given, the file name to save the plot
    """

    alt_err = alt_err[alt_err!=-999]
    if arcsec:
        alt_err *= 3600

    # extract the mean and std
    mu, sigma = norm.fit(alt_err)

    # plot
    fig, ax = plt.subplots(figsize=[15, 10])
    if arcsec:
        ax.set(xlabel=r'Altitude Error ["]', ylabel='Probability', 
               title=r'Tetra3 Altitude Error: $\mu={0:.3f}$", $\sigma={1:.3f}$"'.format(mu, sigma))
    else:
        ax.set(xlabel=r'Altitude Error [$\degree$]', ylabel='Probability', 
               title=r"$\mathrm{Tetra3\ Altitude\ Error:}\ \mu=%.3e\degree,\ \sigma=%.3e\degree$" %(mu, sigma))

    # the histogram of the data
    n, bins, patches = plt.hist(alt_err, bins, density=True, facecolor='darkgreen', 
                                alpha=0.65, rwidth=0.8, edgecolor='black')

    # add a 'best fit' line
    y = norm.pdf(bins, mu, sigma)
    l = plt.plot(bins, y, c='black', ls='--', linewidth=2, label='Gaussian Fit')
    
    plt.legend(loc='upper right')
    plt.grid(True, axis='y')

    if save:
        plt.savefig(save)
    
    #plot
    # plt.show()

    return mu, sigma

def plot_ct(star_ct_col, solve_bool_col, normalize=True, line=False, save=None):
    """
    star_ct_col=the data column of star counts, solve_bool_col=the data column
    describing the solve status, normalize=if True, convert hist data to percentages,
    line=if True, plot lines instead of bar graph, save=if given, the file name
    for the saved plot
    """

    # initialize count holders
    ct_solved_dict = {}
    ct_failed_dict = {}
    ct_timeout_dict = {}

    # fill in the counts
    for ct, sol in zip(star_ct_col, solve_bool_col):
        # ignore low counts
        if ct < 4:
            continue
        
        # initialize ct entry if not defined
        try:
            ct_solved_dict[ct] += 0
        except KeyError:
            ct_solved_dict[ct] = 0
            ct_failed_dict[ct] = 0
            ct_timeout_dict[ct] = 0
            pass
        
        # fill in entries
        if sol == 1:
            ct_solved_dict[ct] += 1
        elif sol == 0:
            ct_failed_dict[ct] += 1
        elif sol == -1:
            ct_timeout_dict[ct] += 1

    # convert dictionaries to count sorted arrays
    ct_solved = []
    for key, value in ct_solved_dict.items():
        ct_solved.append([key, value])
    ct_solved = np.array(ct_solved).T
    _x, _y = ct_solved[0], ct_solved[1]
    ct_solved = np.array([[x, y] for x,y in sorted(zip(_x,_y))], dtype=float).T
    
    ct_failed = []
    for key, value in ct_failed_dict.items():
        ct_failed.append([key, value])
    ct_failed = np.array(ct_failed).T
    _x, _y = ct_failed[0], ct_failed[1]
    ct_failed = np.array([[x, y] for x,y in sorted(zip(_x,_y))], dtype=float).T
    
    ct_timeout = []
    for key, value in ct_timeout_dict.items():
        ct_timeout.append([key, value])
    ct_timeout = np.array(ct_timeout).T
    _x, _y = ct_timeout[0], ct_timeout[1]
    ct_timeout = np.array([[x, y] for x,y in sorted(zip(_x,_y))], dtype=float).T
    
    cts = ct_solved[0]

    if normalize:
        for ct in range(len(cts)):
            try:
                ct_sum = np.sum([ct_solved[1][ct], ct_failed[1][ct], ct_timeout[1][ct]])
                ct_solved[1][ct] = ct_solved[1][ct] / ct_sum
                ct_failed[1][ct] = ct_failed[1][ct] / ct_sum
                ct_timeout[1][ct] = ct_timeout[1][ct] / ct_sum
            except ZeroDivisionError:
                continue

    # plot graph
    fig, ax = plt.subplots(figsize=[20, 10])

    ax.set(xlabel='Star Count', ylabel='{0}Amount Solved/Failed/Timeout'.format('Relative ' if normalize else ''), 
           title='tetra3 Performance for Different Star Counts')
    
    if line:
        plt.plot(*ct_solved, label='Solved', color='green', marker='o', alpha=0.8)
        plt.plot(*ct_failed, label='Failed', color='red', marker='o', alpha=0.75)
        plt.plot(*ct_timeout, label='Timed Out', color='gray', marker='o', alpha=0.8)
    else:
        bar_width = np.min(np.diff(cts))/4
        ct_solved[0] += -bar_width
        plt.bar(*ct_solved, bar_width, label='Solved', color='green', linewidth=2, edgecolor='black', alpha=0.8)
        plt.bar(*ct_failed, bar_width, label='Failed', color='red', linewidth=2, edgecolor='black', alpha=0.75)
        ct_timeout[0] += bar_width
        plt.bar(*ct_timeout, bar_width, label='Timed Out', color='gray', linewidth=2, edgecolor='black', alpha=0.8)
    
    plt.legend(loc='upper right')
    plt.grid(True, axis='y')

    if save:
        plt.savefig(save)

    # plt.show()
    
def analyze_solve_time(star_ct_col, solve_time_col, solve_bool_col, save=None):
    """
    star_ct_col=the data column of star counts, solve_time_col=the data column of
    solve times, solve_bool_col=the data column describing the solve status, 
    normalize=if True, convert hist data to percentages, save=if given, the file 
    name for the saved plot
    """

    # extract relevant data
    star_ct0 = star_ct_col[solve_bool_col == 1]
    solve_time0 = solve_time_col[solve_bool_col == 1]

    star_ct  = np.unique(star_ct0)
    solve_time = np.array([np.mean(solve_time0[star_ct0==ct]) for ct in star_ct])
    solve_time_err = np.array([np.std(solve_time0[star_ct0==ct]) for ct in star_ct])
    
    data = np.array([[x, y] for x,y in sorted(zip(star_ct, solve_time))], dtype=float).T

    # plot
    fig, ax = plt.subplots(figsize=[15, 10])

    ax.set(xlabel='Star Count', ylabel='Solve Time [ms]', 
           title='tetra3 Solve Time vs Star Counts')

    plt.grid()
    
    plt.errorbar(*data, yerr=solve_time_err, fmt='o', capsize=5, c='blue')
    
    if save:
        plt.savefig(save)
    
    # plt.show()

# command line arguments
parser = argparse.ArgumentParser(
    description='Plots the run_sim.py data. Use the options to decide which plots to save. Options include the histogram of the altitude error, the plot of the solve/fail/timeout performance of tetra3 across star counts, and the solve time across star counts.')
parser.add_argument('DATA_NAME', type=str,
                    help='str; the name of the data file, .txt extension assumed')
parser.add_argument('SAVE_NAME', type=str,
                    help='str; the common name of the saved pngs')
parser.add_argument('-e', '--error', action='store_true',
                    help='flag; plot the altitude error histogram')
parser.add_argument('-p', '--perf', action='store_true',
                    help='flag; plot the solve/fail/timeout performance')
parser.add_argument('-t', '--time', action='store_true',
                    help='flag; plot the solve time vs star count')
parser.add_argument('-a', '--all', action='store_true',
                    help='flag; plot all')
args = parser.parse_args()

data_dir = op.join(op.dirname(op.dirname(op.abspath(__file__))), *("data", "images"))
FILE = "{0}/{1}.txt".format(data_dir, args.DATA_NAME)
SAVE = "{0}/{1}".format(data_dir, args.SAVE_NAME)

# obtain header information to print
header = ""
with open(FILE, 'rt') as file:
    for line in file:
        if '#' not in line:
            break
        header += line[2:]
        _temp = line.split()
        if "Total Number" in line: 
            total_number = _temp[5]
        elif "Solved Trials" in line:
            solved_tot = float(_temp[3][:-1])
            solved_per_solve = float(_temp[4][:-1])
            solved_per_tot = float(_temp[8][:-1])
        elif "Solver Fails" in line:
            fail_tot = float(_temp[3][:-1])
            fail_per_solve = float(_temp[4][:-1])
            fail_per_tot = float(_temp[8][:-1])
        elif "Timeout Fails" in line:
            timeout_tot = float(_temp[3][:-1])
            timeout_per_solve = float(_temp[4][:-1])
            timeout_per_tot = float(_temp[8][:-1])
        elif "Low Star Counts" in line:
            low_ct_tot = float(_temp[5][:-1])
            low_ct_per = float(_temp[6][:-1])
        elif "Time" in line:
            tot_time = float(_temp[5])
        elif "skip_header" in line:
            skip_header = int(_temp[3])
        
data = np.genfromtxt(FILE, unpack=True, skip_header=skip_header)

alt_err = data[2]
star_ct = data[4]
solve_bool = data[1]
solve_time = data[9]

print(header)

# plot desired
if args.error or args.all:
    esave = "{0}_alterr.png".format(SAVE)
    analyze_alt_error(alt_err, bins='auto', arcsec=True, save=esave)
if args.perf or args.all:
    psave = "{0}_perf.png".format(SAVE)
    plot_ct(star_ct, solve_bool, normalize=True, line=True, save=psave)
if args.time or args.all:
    tsave = "{0}_time.png".format(SAVE)
    analyze_solve_time(star_ct, solve_time, solve_bool, save=tsave)
