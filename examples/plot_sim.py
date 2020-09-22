"""
Plots the data obtained from run_sim.py or run_sim.sh
"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams.update({'font.size': 20})
import argparse
import glob
import os.path as op
data_dir = op.join(op.dirname(op.dirname(op.abspath(__file__))), *("data", "sim"))

def get_header(filename):
    """Return the header and header dict of run_sim.py data logs"""
    header = ""
    header_dict = {}
    with open(filename, 'rt') as file:
        for line in file:
            if '#' not in line:
                break
            header += line[2:]
            _temp = line.split()
            if "Total Number" in line:
                header_dict['total_number'] = _temp[5]
            elif "Solved Trials" in line:
                header_dict['solved_tot'] = float(_temp[3][:-1])
                header_dict['solved_per_solve'] = float(_temp[4][:-1])
                header_dict['solved_per_tot'] = float(_temp[8][:-1])
            elif "Solver Fails" in line:
                header_dict['fail_tot'] = float(_temp[3][:-1])
                header_dict['fail_per_solve'] = float(_temp[4][:-1])
                header_dict['fail_per_tot'] = float(_temp[8][:-1])
            elif "Timeout Fails" in line:
                header_dict['timeout_tot'] = float(_temp[3][:-1])
                header_dict['timeout_per_solve'] = float(_temp[4][:-1])
                header_dict['timeout_per_tot'] = float(_temp[8][:-1])
            elif "Low Star Counts" in line:
                header_dict['low_ct_tot'] = float(_temp[5][:-1])
                header_dict['low_ct_per'] = float(_temp[6][:-1])
            elif "Time" in line:
                header_dict['tot_time'] = float(_temp[5])
            elif "skip_header" in line:
                header_dict['skip_header'] = int(float(_temp[3]))
            elif "     " in line:
                header_dict['columns'] = _temp[1:]
    return header, header_dict

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
        plt.savefig(save, bbox_inches='tight', pad_inches=0)
    
    #plot
    # plt.show()

    return mu, sigma

def plot_ct(star_ct_col, solve_bool_col, normalize=True, line=False, save=None, plot=True):
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

    if not plot:
        return cts, ct_solved[1], ct_failed[1], ct_timeout[1]
    else:
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
            plt.savefig(save, bbox_inches='tight', pad_inches=0)

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
        plt.savefig(save, bbox_inches='tight', pad_inches=0)
    
    # plt.show()

def analyze_param(param: str, filenames, savename, p_alt_err: bool = True, p_perf: bool = True, p_time: bool = True):
    """Plot Alt_err vs param, Perf vs star ct vs param, and solve time vs star ct vs param"""

    if param=='rotation':
        param_unit = ' ["/s]'
    else:
        param_unit = "" 
    data_dict = {}
    rates = {}
    par_vals = []
    for FILE in filenames:
         hdr_dict = get_header(FILE)[1]
         solve_rate = hdr_dict['solved_per_solve']
         fail_rate = hdr_dict['fail_per_solve']
         timeout_rate = hdr_dict['timeout_per_solve']
         
         data = np.genfromtxt(FILE, unpack=True, skip_header=hdr_dict['skip_header'])
         alt_err = data[2]
         star_ct = data[4]
         solve_bool = data[1]
         solve_time = data[9]

         param_ind = hdr_dict['columns'].index(param)
         par_val = data[param_ind]
         par_val = par_val[par_val!=-999][0]

         data_dict[str(par_val)] = np.array([alt_err, star_ct, solve_bool, solve_time])
         rates[str(par_val)] = np.array([solve_rate, fail_rate, timeout_rate])
         par_vals.append(par_val)

    par_vals = np.array(par_vals)

    def param_alt_err(data_dict, par_vals, save):
        """Plot alt_err vs param"""

        # extract the mean and std
        altdata = []
        for par in data_dict:
            # append [mu, sigma]
            alt = data_dict[par][0]
            alt = alt[alt!=-999]*60*60
            altdata.append(np.array(norm.fit(alt)))
        altdata = np.array(altdata).T    
            
        # plot
        fig, ax = plt.subplots(figsize=[15, 10])

        ax.set(xlabel=param+param_unit, ylabel='Altitude Error ["]', 
               title='Altitude Error vs {0}'.format(param))
        if param=='snr':
            plt.xscale('log')

        if param=='rotation':
            mean_rot = 152
            plt.axvline(x=mean_rot, c='red', alpha=0.8, label='ANITA Mean Rotation')
            plt.legend(loc='upper right')

        plt.errorbar(par_vals, altdata[0], yerr=altdata[1], fmt='o', capsize=5, c='black', alpha=0.8)
        plt.grid(True)

        if save:
            plt.savefig(save, bbox_inches='tight', pad_inches=0.1)

        #plot
        # plt.show()

    def param_perf(data_dict, rates, par_vals, savename):

        # extract rates
        ratedata = []
        for par in rates:
            ratedata.append(np.array(rates[par]))
        ratedata = np.array(ratedata).T
        plotdata = np.array([[p, s, f, t] for p,s,f,t in sorted(zip(par_vals,*ratedata))], dtype=float).T
        
        # plot rates
        fig, ax = plt.subplots(figsize=[15, 10])

        ax.set(xlabel=param+param_unit, ylabel='Percentage Solved/Failed/Timeout', 
               title='Solve Rates vs {0}'.format(param))
        if param=='snr':
            plt.xscale('log')

        x = plotdata[0]
        plt.plot(x, plotdata[1], label='Solved', color='green', marker='o', alpha=0.8)
        plt.plot(x, plotdata[2], label='Failed', color='red', marker='o', alpha=0.75)
        plt.plot(x, plotdata[3], label='Timed Out', color='gray', marker='o', alpha=0.8)

        plt.grid(True)
        plt.legend(loc='upper right')

        if savename:
            save = "{0}_perfrate.png".format(savename)
            plt.savefig(save, bbox_inches='tight', pad_inches=0.1)

        #plot
        # plt.show()

        # plot solve performance
        fig2, ax2 = plt.subplots(figsize=[15, 10])
        ax2.set(xlabel='Star Count', ylabel='Relative Amount Solved', 
                title='Solve Performance vs Star Count at Different {0}'.format(param))
        plt.grid(True)

        # set color map
        evenly_spaced_interval = np.linspace(0, 1, len(x))
        colors = [cm.rainbow(i) for i in evenly_spaced_interval]

        for par in range(len(x)):
            ind = str(float(x[par]))
            star_ct_col = data_dict[ind][1]
            solve_bool_col = data_dict[ind][2]
            xdata, ydata, _, _ = plot_ct(star_ct_col, solve_bool_col, normalize=True, plot=False)
            plt.plot(xdata, ydata, label=ind, marker='o', color=colors[par], alpha=0.8)
        handles, labels = ax2.get_legend_handles_labels()
        plt.legend(handles[::-1], labels[::-1], title=param+param_unit, loc='lower right')
        
        if savename:
            save = "{0}_perf.png".format(savename)
            plt.savefig(save, bbox_inches='tight', pad_inches=0.1)

    def param_time(data_dict, par_vals, savename):

        fig, ax = plt.subplots(figsize=[15, 10])
        ax.set(xlabel='Star Count', ylabel='Solve Time [ms]', 
               title='Solve Time vs Star Counts at Different {0}'.format(param))
        plt.grid()

        # set color map
        evenly_spaced_interval = np.linspace(0, 1, len(par_vals))
        colors = [cm.rainbow(i) for i in evenly_spaced_interval]

        par_vals = sorted(par_vals)
        
        for par in (0, -1):

            ind = str(float(par_vals[par]))
            star_ct_col = data_dict[ind][1]
            solve_time_col = data_dict[ind][3]
            solve_bool_col = data_dict[ind][2]
            
            # extract relevant data
            star_ct0 = star_ct_col[solve_bool_col == 1]
            solve_time0 = solve_time_col[solve_bool_col == 1]

            star_ct  = np.unique(star_ct0)
            solve_time = np.array([np.mean(solve_time0[star_ct0==ct]) for ct in star_ct])
            solve_time_err = np.array([np.std(solve_time0[star_ct0==ct]) for ct in star_ct])

            data = np.array([[x, y] for x,y in sorted(zip(star_ct, solve_time))], dtype=float).T

            plt.errorbar(*data, yerr=solve_time_err, fmt='o', capsize=5, label=ind, c=colors[par])

        handles, labels = ax.get_legend_handles_labels()
        plt.legend(handles[::-1], labels[::-1], title=param+param_unit, loc='upper right')
            
        if savename:
            plt.savefig(savename, bbox_inches='tight', pad_inches=0.1)

            
    if p_alt_err:
        SAVE = "{0}_alterr.png".format(savename)
        param_alt_err(data_dict, par_vals, SAVE)

    if p_perf:
        param_perf(data_dict, rates, par_vals, savename)

    if p_time:
        SAVE = "{0}_time.png".format(savename)
        param_time(data_dict, par_vals, SAVE)
    
# command line arguments
parser = argparse.ArgumentParser(
    description='Plots the run_sim.py data. Use the options to decide which plots to save or none to plot all. Options include the histogram of the altitude error, the plot of the solve/fail/timeout performance of tetra3 across star counts, and the solve time across star counts. Use --param to perform the analysis across multiple data with a varying specific parameter whose filenames contain a common string (ie. SNR10.txt, SNR1.txt, SNR0.1.txt,... -> pass SNR* as DATA_NAME). Note that the data is assumed to be in the data/sim/ directory.')
parser.add_argument('DATA_NAME', type=str,
                    help='str; the name of the data file, .txt extension assumed. If --param, this is the file name pattern e.g. SNR*, ROT*ANG*, etc.')
parser.add_argument('SAVE_NAME', type=str,
                    help='str; the common name of the saved pngs')
parser.add_argument('-e', '--error', action='store_true',
                    help='flag; plot the altitude error histogram')
parser.add_argument('-p', '--perf', action='store_true',
                    help='flag; plot the solve/fail/timeout performance')
parser.add_argument('-t', '--time', action='store_true',
                    help='flag; plot the solve time vs star count')
parser.add_argument('--param', type=str, choices=['none','snr','rotation'], default='none',
                    help='str; analyze list of data files differing in this specific parameter')
args = parser.parse_args()

if args.param == 'none':

    plot_all = False if args.error or args.perf or args.time else True
    
    FILE = "{0}/{1}.txt".format(data_dir, args.DATA_NAME)
    SAVE = "{0}/{1}".format(data_dir, args.SAVE_NAME)

    hdr, hdr_dict = get_header(FILE)
    
    data = np.genfromtxt(FILE, unpack=True, skip_header=hdr_dict['skip_header'])

    alt_err = data[2]
    star_ct = data[4]
    solve_bool = data[1]
    solve_time = data[9]

    print(hdr)

    # plot desired
    if args.error or plot_all:
        esave = "{0}_alterr.png".format(SAVE)
        analyze_alt_error(alt_err, bins='auto', arcsec=True, save=esave)
    if args.perf or plot_all:
        psave = "{0}_perf.png".format(SAVE)
        plot_ct(star_ct, solve_bool, normalize=True, line=True, save=psave)
    if args.time or plot_all:
        tsave = "{0}_time.png".format(SAVE)
        analyze_solve_time(star_ct, solve_time, solve_bool, save=tsave)

else:

    filenames = glob.glob("{0}/{1}".format(data_dir, args.DATA_NAME))
    print(filenames)
    SAVE = "{0}/{1}.txt".format(data_dir, args.SAVE_NAME)

    # p_alt_err = False if (args.error or args.perf) and not args.error else True
    # p_perf = False if (args.error or args.perf) and not args.perf else True
    # p_time = args.time
    
    # analyze_param(args.param, filenames, SAVE, p_alt_err=p_alt_err, p_perf=p_perf, p_time=p_time)
