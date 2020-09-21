"""
This script combines a number of result files with a
common basename (RESULT1.txt, RESULT2.txt, ...) into
one result file of a given name.
"""

import os.path as op
data_dir = op.join(op.dirname(op.dirname(op.abspath(__file__))), *("data", "sim"))
import argparse

# command line arguments
parser = argparse.ArgumentParser(
    description='Combines Monte Carlo Sim data files from run_sim.py into one.')
parser.add_argument('TRIAL_NAME', type=str,
                    help='str; data base name')
parser.add_argument('COMBINE_NAME', type=str,
                    help='str; name of combined file')
parser.add_argument('-s', '--START', type=int, default=1,
                    help='int; data file start number')
parser.add_argument('-N', '--N_TRIALS', type=int, default=10,
                    help='int; number of runs')

args = parser.parse_args()

TN = args.TRIAL_NAME
CN = args.COMBINE_NAME
START = args.START
END = args.N_TRIALS + START

# Initialize combined header information
header = ''
total_number = 0.
solved_tot = 0.
solved_per_solve = 0.
solved_per_tot = 0.
fail_tot = 0.
fail_per_solve = 0.
fail_per_tot = 0.
timeout_tot = 0.
timeout_per_solve = 0.
timeout_per_tot = 0.
low_ct_tot = 0.
low_ct_per = 0.
tot_time = 1000000.
snr = None
skip = 0.

# extract header data from each file
for i in range(START, END):
    FILE = '{0}/{1}{2}.txt'.format(data_dir,TN,i)
    try:
        with open(FILE, 'rt') as file:
            for line in file:
                if '#' not in line:
                    break
                header += line[2:]
                _temp = line.split()
                if "Total Number" in line: 
                    _total_number = float(_temp[5])
                elif "Solved Trials" in line:
                    _solved_tot = float(_temp[3][:-1])
                elif "Solver Fails" in line:
                    _fail_tot = float(_temp[3][:-1])
                elif "Timeout Fails" in line:
                    _timeout_tot = float(_temp[3][:-1])
                elif "Low Star Counts" in line:
                    _low_ct_tot = float(_temp[5][:-1])
                elif "Time" in line:
                    _tot_time = float(_temp[5])
                elif "SNR:" in line:
                    snr = float(_temp[2])
                elif "skip_header" in line:
                    skip = int(float(_temp[3]))

            # update combined header information
            total_number += _total_number
            solved_tot += _solved_tot
            fail_tot += _fail_tot
            timeout_tot += _timeout_tot
            low_ct_tot += _low_ct_tot
            tot_time = min(tot_time, _tot_time)
            
    except FileNotFoundError:
        continue

# more combined header information
solver_tries = total_number - low_ct_tot
solved_per_solve = 0 if solver_tries == 0 else 100 * solved_tot / solver_tries
solved_per_tot = 0 if total_number == 0 else 100 * solved_tot / total_number
fail_per_solve = 0 if solver_tries == 0 else 100 * fail_tot / solver_tries
fail_per_tot = 0 if total_number == 0 else 100 * fail_tot / total_number
timeout_per_solve = 0 if solver_tries == 0 else 100 * timeout_tot / solver_tries
timeout_per_tot = 0 if total_number == 0 else 100 * timeout_tot / total_number
low_ct_per = 0 if total_number == 0 else 100 * low_ct_tot / total_number
snr_str = "# SNR: {0}\n".format(snr) if snr else ""

OUT = '{0}/{1}.txt'.format(data_dir,CN)
with open(OUT, 'w') as out:

    # write the header
    out.write(
            '# STEREO Data for run titled "{0}"\n'.format(CN) + 
            '#\n# Total Number of Trials: {0}\n'.format(total_number) +
            '# Solved Trials: {0}, {1}% (of attempted solves), {2}% (total)\n'.format(
                solved_tot, solved_per_solve, solved_per_tot) +
            '# Solver Fails: {0}, {1}% (of attempted solves), {2}% (total)\n'.format(
                fail_tot, fail_per_solve, fail_per_tot) +
            '# Timeout Fails: {0}, {1}% (of attempted solves), {2}% (total)\n'.format(
                timeout_tot, timeout_per_solve, timeout_per_tot) +
            '# Low Star Counts (<{0}): {1}, {2}%\n'.format(
                4, low_ct_tot, low_ct_per) +
            '# Total Data Collection Time: {0} hrs\n'.format(tot_time) +
            snr_str +
            '# skip_header = {0}\n'.format(skip) +
            '#\n# Column Notes\n' + 
            '# num: Data Index\n' + 
            '# solved: 1 = Success, 0 = Fail, -1 = Undetermined (low star count or timed out)\n' + 
            '# alt_error: Error in Altitude (from 45 deg)\n' + 
            '# A4_index: Index used for ANITA4 data\n' + 
            '# star_ct: Number of Stars in the field\n' + 
            '# real_ra: Real RA of the center of the field in degrees\n' +
            '# real_dec: Real DEC of the center of the field in degrees\n' +
            '# ra: Result RA of the center of the field in degrees\n' + 
            '# dec: Result DEC of the center of the field in degrees\n' +
            '# solve_time: Total time taken to solve (or fail) in seconds\n' +
            '#\n# num     solved     alt_error     A4_index     star_ct     real_ra     real_dec     ra     dec     solve_time\n'
        )

    # write the data from each file
    n = 1
    for i in range(START, END):
        FILE = '{0}/{1}{2}.txt'.format(data_dir,TN,i)
        try:
            with open(FILE, 'rt') as file:
                for line in file:
                    if '#' not in line:
                        newline = str(n) + line[line.find(' '):]
                        out.write(newline)
                        n += 1
        except FileNotFoundError:
            continue
