"""
This script runs the Monte Carlo Sim:
1. A flightpath from ANITA is chosen
2. An image is generated using the ANITA data and assumes a
   center altitude coordinate of 45 degrees (random azimuth)
3. The image is run through the tetra3 star tracking algorithm
4. The solution's center coords are compared to the known

The data is written in a txt file of a given name. Use
python run_sim.py -h
to see other parameters to change on the command line.

For a parallel process of this simulation, use run_sim.sh .
"""

import stereo
from stereo.sim import *
from stereo.hardware import ASI183MM_Sigma135mm as cam
import stereo.flightpath as flightpath

import tetra3

import argparse
import itertools
import logging
import os.path as op
data_dir = op.join(op.dirname(op.dirname(op.abspath(__file__))), *("data", "sim"))
import random
from time import perf_counter as timestamp


import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astropy.utils.iers.iers import conf
conf.iers_auto_url = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
conf.auto_download = False


def run_sim(
        N: int,
        run_name: str,
        output: str = None,
        star_min: int = 4,
        snr: float = None,
        rotation: float = None,
        verbose_factor: int = 100,
        single_run: bool = False,
        save_image: bool = False,
        plot: bool = False,
        exit_fail: int = 100,
) -> None:
    """
    Perform STEREO Monte Carlo.

    Parameters
    ----------
    N : int
        Desired number of successfully solved trials.
    run_name : str
        A name for this run.
    output : str
        The file name of the data log. If not given, no log is saved.
    star_min : int
        Skips the trial if the generated field has less than this many stars.
    snr : float
        If given, forces the image to the desired signal-to-noise ratio.
    rotation : float
        The rotation in arcsec/s. If None, rotation is selected randomly from
        a Gaussian range in ANITA flightpath data.
    verbose_factor : int
        Prints the progress every this many trials.
    single_run : bool
        If True, ignores N and runs for exactly one trial, successful or not.
        Prints the results immediately.
    save_image : bool
        If True, save the image (as a compressed fits image).
    plot : bool
        If True, immediately plots the image.
    exit_fail : int
        How many images in a row to fail to solve or timeout before giving up.
    """

    # sim intialization
    data_pts = 0
    trial_num = 0
    star_num_fail = 0
    solver_fail = 0
    timeout_fail = 0
    exit_fail_ct = 0
    data_generator = ()
    
    # image and solver params
    cen_alt = 45 #degrees
    mag_limit = 8
    psf_sigma = 3
    centroid_dict = {
        'sigma':3, 
        'image_th':None, 
        'filtsize':10, 
        'binary_open':True, 
        'bg_sub_mode':'global_median', 
        'sigma_mode':'global_median_abs', 
        'centroid_window':None, 
        'max_area':None, 
        'min_area':5, 
        'max_sum':None, 
        'min_sum':100, 
        'max_axis_ratio':None,
        'downsample':4
    }
    fov = cam.fov
    solve_dict = {
        'fov_estimate':fov[0], 
        'fov_max_error':0.1, 
        'pattern_checking_stars':6,
        'match_radius':0.01,
        'match_threshold':1e-8
    }
    anita = flightpath.load_flight(4)
    db = "db_fov6ps15cs20pme1"
    t3 = tetra3.Tetra3()
    t3._logger.setLevel(logging.ERROR)
    t3.load_database(db)
    t3.database_properties['pattern_max_error'] = 0.005

    img_dir = op.join(op.dirname(data_dir), "images")

    if not single_run:
        print("{0}: Beginning Data Collection...".format(run_name))
    
    total_t0 = timestamp()
    
    try:

        # run until N non-None alt_err data points
        while data_pts < N:

            # begin new trial
            trial_num += 1
            star_ct = 0

            # UI progress display
            if int(verbose_factor) == 0 or trial_num % int(verbose_factor) == 0:
                print("\r{0}: Solving Trial ".format(run_name), trial_num, "... \n(", data_pts,
                      " Solved Successfully, ", star_num_fail, " Star Count Fails, ", solver_fail,
                      " Failed to Solve, ", timeout_fail, " Timed Out)", sep="")

            # handle repeated failures
            if exit_fail_ct >= exit_fail:
                print("{0}: Failed {1} times in a row. Stopping sim...".format(run_name, exit_fail))
                raise KeyboardInterrupt

            # get a random center coord from ANITA flight
            ind = random.randrange(0, 37279)
            lat = anita.latitude[ind]
            lon = anita.longitude[ind]
            alt = anita.altitude[ind]
            time = Time(anita.realTime[ind], format="unix")
            obs = Observation(lat=lat, lon=lon, alt=alt, time=time)
            altaz = obs.altaz_frame
            cen_az = random.randrange(0, 360) # degrees
            center_radec = SkyCoord(alt=cen_alt, az=cen_az, unit=u.deg, frame=altaz).transform_to('icrs')

            # generate the image
            image = generate_starfield(center_radec, cam, obs, mag_limit, psf_sigma)

            # check number of stars
            star_ct = len(image.image_stars)

            # star count test
            if star_ct < star_min:
                star_num_fail += 1
                
                # add data line
                trials = trial_num
                solved = -1
                alt_err = -999
                ind_val = ind
                num_stars = star_ct
                real_ra = center_radec.ra.value
                real_dec = center_radec.dec.value
                result_ra = -999
                result_dec = -999
                solve_t = -999
                data_snr = -999
                data_rot = -999
            
                data_line = [[
                    trials, solved, alt_err, ind_val, num_stars, 
                    real_ra, real_dec, result_ra, result_dec, solve_t,
                    data_snr, data_rot
                ]]
                data_generator = itertools.chain(data_generator, data_line)
                
                continue

            # save image
            if save_image:
                img = fits.CompImageHDU(image.image)
                filename = '{0}/{1}_{2:04}.fits'.format(img_dir, run_name, trial_num)
                img.writeto(filename, overwrite=True)

            # make the image
            image.rotation = rotation if rotation else np.random.normal(loc=152.504, scale=949.721)
            image.snr = snr
            make_image = image.image

            # plot image
            if plot:
                image.plot()
                
            # run star tracker algorithm
            try:
                t3result = run_tetra3(image, database=db, t3=t3, return_result=True,
                                      **solve_dict, **centroid_dict)
                solve_tf = t3result['T_solve']
                if t3result["RA"] is None or t3result["Dec"] is None:
                    result = None
                else:
                    result = SkyCoord(ra=t3result["RA"], dec=t3result["Dec"], unit=u.deg)
            except StopIteration:
                timeout_fail += 1
                exit_fail_ct += 1
                
                # add data line
                trials = trial_num
                solved = -1
                alt_err = -999
                ind_val = ind
                num_stars = star_ct
                real_ra = center_radec.ra.value
                real_dec = center_radec.dec.value
                result_ra = -999
                result_dec = -999
                solve_t = -999
                data_snr = image.snr
                data_rot = image.rotation
            
                data_line = [[
                    trials, solved, alt_err, ind_val, num_stars, 
                    real_ra, real_dec, result_ra, result_dec, solve_t,
                    data_snr, data_rot
                ]]
                data_generator = itertools.chain(data_generator, data_line)

                if single_run:
                    break
                else:
                    continue
            except KeyboardInterrupt:
                raise
            except:
                print("Mystery error occurred at")
                print('Index: {0}, Center: {1}'.format(ind, center_radec))
                continue
            
            if not result: # no solution
                solver_fail += 1
                exit_fail_ct += 1
                
                # add data line
                trials = trial_num
                solved = 0
                alt_err = -999
                ind_val = ind
                num_stars = star_ct
                real_ra = center_radec.ra.value
                real_dec = center_radec.dec.value
                result_ra = -999
                result_dec = -999
                solve_t = solve_tf
                data_snr = image.snr
                data_rot = image.rotation
            
                data_line = [[
                    trials, solved, alt_err, ind_val, num_stars, 
                    real_ra, real_dec, result_ra, result_dec, solve_t,
                    data_snr, data_rot
                ]]
                data_generator = itertools.chain(data_generator, data_line)

                if single_run:
                    break
                else:
                    continue

            else: # success
                data_pts += 1
                exit_fail_ct = 0 #reset
                
                # get alt error
                result_altaz = result.transform_to(altaz)
                result_alt = result_altaz.alt.value
                err = cen_alt - result_alt
                
                # add data line
                trials = trial_num
                solved = 1
                alt_err = err
                ind_val = ind
                num_stars = star_ct
                real_ra = center_radec.ra.value
                real_dec = center_radec.dec.value
                result_ra = result.ra.value
                result_dec = result.dec.value
                solve_t = solve_tf
                data_snr = image.snr
                data_rot = image.rotation
            
                data_line = [[
                    trials, solved, alt_err, ind_val, num_stars, 
                    real_ra, real_dec, result_ra, result_dec, solve_t,
                    data_snr, data_rot
                ]]
                data_generator = itertools.chain(data_generator, data_line)

                if single_run:
                    break
                else:
                    continue
                
        total_t = (timestamp() - total_t0) / 60 / 60

        if single_run:
            print('RA: {7}\nDEC: {8}\nTime: {9}'.format(*next(data_generator)))

        if output:
            
            outname = "{0}/{1}.txt".format(data_dir,output)
            out = open(outname, 'w')
            
            print("{0}: Done. Saving data to {1}...".format(run_name, output))

            # header data
            data_pts_per_solver = 100 * data_pts / (trial_num - star_num_fail)
            data_pts_per_total = 100 * data_pts / trial_num
            solver_fail_per_solver = 100 * solver_fail / (trial_num - star_num_fail)
            solver_fail_per_total = 100 * solver_fail / trial_num
            timeout_fail_per_solver = 100 * timeout_fail / (trial_num - star_num_fail)
            timeout_fail_per_total = 100 * timeout_fail / trial_num
            star_num_fail_per_total = 100 * star_num_fail / trial_num
            skip_str = '# skip_header = 23\n'

            out.write(
                '# STEREO Data for run titled "{0}"\n'.format(run_name) + 
                '#\n# Total Number of Trials: {0}\n'.format(trial_num) +
                '# Solved Trials: {0}, {1}% (of attempted solves), {2}% (total)\n'.format(
                    data_pts, data_pts_per_solver, data_pts_per_total) +
                '# Solver Fails: {0}, {1}% (of attempted solves), {2}% (total)\n'.format(
                    solver_fail, solver_fail_per_solver, solver_fail_per_total) +
                '# Timeout Fails: {0}, {1}% (of attempted solves), {2}% (total)\n'.format(
                    timeout_fail, timeout_fail_per_solver, timeout_fail_per_total) +
                '# Low Star Counts (<{0}): {1}, {2}%\n'.format(
                    star_min, star_num_fail, star_num_fail_per_total) +
                '# Total Data Collection Time: {0} hrs\n'.format(total_t) + 
                skip_str + 
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
                '# snr: mean(signal) / sqrt(mean(background))\n' +
                '# rotation: payload/camera rotation in arcsec/s\n' +
                '#\n# num     solved     alt_error     A4_index     star_ct     real_ra     real_dec     ra     dec     solve_time     snr     rotation\n'
            )

            for data in data_generator:
                out.write(
                    '{0}     {1}     {2}     {3}     {4}     {5}     {6}     {7}     {8}     {9}     {10}     {11}\n'.format(*data)
                )

            out.close()
        
        else:
            print("{0}: Done.".format(run_name))
            
        # End
    
    # return what's there so far if interrupted
    except:

        if output:

            outname = "{0}/{1}.txt".format(data_dir,output)
            out = open(outname, 'w')

            total_t = (timestamp() - total_t0) / 60 / 60

            print("{0}: Data collection interrupted. Saving data to {1}...".format(run_name, output))

            # header data
            data_pts_per_solver = 100 * data_pts / (trial_num - star_num_fail)
            data_pts_per_total = 100 * data_pts / trial_num
            solver_fail_per_solver = 100 * solver_fail / (trial_num - star_num_fail)
            solver_fail_per_total = 100 * solver_fail / trial_num
            timeout_fail_per_solver = 100 * timeout_fail / (trial_num - star_num_fail)
            timeout_fail_per_total = 100 * timeout_fail / trial_num
            star_num_fail_per_total = 100 * star_num_fail / trial_num
            skip_str = '# skip_header = 23\n'

            out.write(
                '# STEREO Data for run titled "{0}"\n'.format(run_name) + 
                '#\n# Total Number of Trials: {0}\n'.format(trial_num) +
                '# Solved Trials: {0}, {1}% (of attempted solves), {2}% (total)\n'.format(
                    data_pts, data_pts_per_solver, data_pts_per_total) +
                '# Solver Fails: {0}, {1}% (of attempted solves), {2}% (total)\n'.format(
                    solver_fail, solver_fail_per_solver, solver_fail_per_total) +
                '# Timeout Fails: {0}, {1}% (of attempted solves), {2}% (total)\n'.format(
                    timeout_fail, timeout_fail_per_solver, timeout_fail_per_total) +
                '# Low Star Counts (<{0}): {1}, {2}%\n'.format(
                    star_min, star_num_fail, star_num_fail_per_total) +
                '# Total Data Collection Time: {0} hrs\n'.format(total_t) + 
                skip_str + 
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
                '# snr: mean(signal) / sqrt(mean(background))\n' +
                '# rotation: payload/camera rotation in arcsec/s\n' +
                '#\n# num     solved     alt_error     A4_index     star_ct     real_ra     real_dec     ra     dec     solve_time     snr     rotation\n'
            )

            for data in data_generator:
                out.write(
                    '{0}     {1}     {2}     {3}     {4}     {5}     {6}     {7}     {8}     {9}     {10}     {11}\n'.format(*data)
                )

            out.close()

        else:
            print("{0}: Interrupted.".format(run_name))
        
        # End

# command line arguments
parser = argparse.ArgumentParser(
    description='Performs a Monte Carlo of tetra3 over generated images using ANITA data for a given number of solved trials.')
parser.add_argument('RUN_NAME', type=str,
                    help='str; output filename of results')
parser.add_argument('-t', '--TOTAL_SOLVE_TRIALS', type=float, default=100,
                    help='int; number of solved trials to accomplish; default=100')
parser.add_argument('--SNR', type=float, default=None,
                    help='float; the desired signal-to-noise ratio; default=None')
parser.add_argument('--ROTATION', type=float, default=None,
                    help='float; the desired rotation of the payload')
parser.add_argument('-v', '--VERBOSE_FACTOR', type=int, default=100,
                    help='int; print progress per this factor; default=100')
parser.add_argument('-s', '--SINGLE_RUN', action='store_true',
                    help='perform for only a single image')
parser.add_argument('-i', '--SAVE_IMAGE', action='store_true',
                    help='save the image (as a compressed fits image)')
parser.add_argument('-p', '--PLOT', action='store_true',
                    help='immediately plot the image')

args = parser.parse_args()

TST = args.TOTAL_SOLVE_TRIALS
RN = args.RUN_NAME
SNR = args.SNR
ROT = args.ROTATION
VF = args.VERBOSE_FACTOR
SR = args.SINGLE_RUN
SI = args.SAVE_IMAGE
P = args.PLOT
output = None if SR else RN

run_sim(
    TST,
    RN,
    star_min=4,
    snr=SNR,
    rotation=ROT,
    verbose_factor=VF,
    single_run=SR,
    save_image=SI,
    plot=P,
    exit_fail=100,
    output=output
)
