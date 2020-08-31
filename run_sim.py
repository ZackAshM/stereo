import stereo
from stereo.sim import *
from stereo.hardware import ASI183MM_Sigma135mm as cam

import tetra3

import itertools

from time import perf_counter as timestamp

from astropy.io import fits

import numpy as np

import random

from astropy.coordinates import SkyCoord
from astropy import units as u

def run_sim(N, trial_name, output = None, star_min = 4, verbose_factor=100):
    """N = number of non-None points, trial = trial name, output = filename of output,
    star_min = min star count for solve attempt, verbose_factor = print trial progress"""

    # sim intialization
    data_pts = 0
    trial_num = 0
    star_num_fail = 0
    solver_fail = 0
    timeout_fail = 0
    data_generator = ()
    
    # image and solver params
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
        'max_axis_ratio':None
    }
    fov = cam.fov
    fov_err = fov[2] - min(fov[0], fov[1])
    solve_dict = {
        'fov_estimate':fov[0], 
        'fov_max_error':fov_err, 
        'pattern_checking_stars':-1, 
        'match_radius':0.01, 
        'match_threshold':1e-9
    }
    db = "db_fov10ps15cs20pme01"
    ds = 4
    t3 = tetra3.Tetra3()
    t3.load_database(db)
    
    # open output
    if output:
        out = open(output, 'w')

    print("Beginning Data Collection...")
        
    total_t0 = timestamp()
    
    try:

        # run until N non-None alt_err data points
        while data_pts < N:

            # begin new trial
            trial_num += 1
            star_ct = 0

            # UI progress display
            if trial_num % int(verbose_factor) == 0:
                print("Solving Trial ", trial_num, "... \n(", data_pts, " Solved Successfully, ", 
                      star_num_fail, " Star Count Fails, ", solver_fail, " Failed to Solve, ", 
                      timeout_fail, " Timed Out)", sep="")

            # get a random center coord from ANITA flight
            ind = random.randrange(0, 37279)
            obs = load_anita_observation(4, ind)
            cen_alt = 45 #degrees                                                                                                           
            cen_az = random.randrange(0, 360) # degrees                                                                                     
            altaz = obs.altaz_frame
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
            
                data_line = [[
                    trials, solved, alt_err, ind_val, num_stars, 
                    real_ra, real_dec, result_ra, result_dec, solve_t
                ]]
                data_generator = itertools.chain(data_generator, data_line)
                
                continue
                
            if star_ct >= 10:
                img = fits.CompImageHDU(image.image)
                filename = 'data/images/{0}_{1:04}.fits'.format(trial_name, trial_num)
                img.writeto(filename, overwrite=True)

            # run star tracker algorithm
            solve_t0 = timestamp()
            try:
                result = run_tetra3(image, database=db, t3=t3, downsample=ds, **solve_dict, **centroid_dict)
            except StopIteration:
                timeout_fail += 1
                
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
            
                data_line = [[
                    trials, solved, alt_err, ind_val, num_stars, 
                    real_ra, real_dec, result_ra, result_dec, solve_t
                ]]
                data_generator = itertools.chain(data_generator, data_line)
                
                continue
            
            solve_tf = timestamp() - solve_t0

            if not result: # no solution
                solver_fail += 1
                
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
            
                data_line = [[
                    trials, solved, alt_err, ind_val, num_stars, 
                    real_ra, real_dec, result_ra, result_dec, solve_t
                ]]
                data_generator = itertools.chain(data_generator, data_line)
                
                continue

            else: # success
                data_pts += 1
                
                # success, but star count low? Save it
                if star_ct < 10:
                    img = fits.CompImageHDU(image.image)
                    filename = 'data/images/{0}_{1:04}.fits'.format(trial_name, trial_num)
                    img.writeto(filename, overwrite=True)

                # get alt error
                result_altaz = result.transform_to(altaz)
                result_alt = result_altaz.alt.value
                err = np.abs(cen_alt - result_alt)
                
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
            
                data_line = [[
                    trials, solved, alt_err, ind_val, num_stars, 
                    real_ra, real_dec, result_ra, result_dec, solve_t
                ]]
                data_generator = itertools.chain(data_generator, data_line)
                
        total_t = (timestamp() - total_t0) / 60 / 60

        print("Done. Saving data to {0}...".format(output))
        
        # header data
        data_pts_per_solver = 100 * data_pts / (trial_num - star_num_fail)
        data_pts_per_total = 100 * data_pts / trial_num
        solver_fail_per_solver = 100 * solver_fail / (trial_num - star_num_fail)
        solver_fail_per_total = 100 * solver_fail / trial_num
        timeout_fail_per_solver = 100 * timeout_fail / (trial_num - star_num_fail)
        timeout_fail_per_total = 100 * timeout_fail / trial_num
        star_num_fail_per_total = 100 * star_num_fail / trial_num
        
        out.write(
            '# STEREO Data for run titled "{0}"\n'.format(trial_name) + 
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
            '# skip_header = 23\n' + 
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
        
        for data in data_generator:
            out.write(
                '{0}     {1}     {2}     {3}     {4}     {5}     {6}     {7}     {8}     {9}\n'.format(*data)
            )

        out.close()
            
        # End
    
    # return what's there so far if interrupted
    except KeyboardInterrupt:

        total_t = (timestamp() - total_t0) / 60 / 60
        
        print("Data collection interrupted. Saving data to {0}...".format(output))
        
        data_pts_per_solver = 100 * data_pts / (trial_num - star_num_fail)
        data_pts_per_total = 100 * data_pts / trial_num
        solver_fail_per_solver = 100 * solver_fail / (trial_num - star_num_fail)
        solver_fail_per_total = 100 * solver_fail / trial_num
        timeout_fail_per_solver = 100 * timeout_fail / (trial_num - star_num_fail)
        timeout_fail_per_total = 100 * timeout_fail / trial_num
        star_num_fail_per_total = 100 * star_num_fail / trial_num
        
        out.write(
            '# STEREO Data for run titled "{0}"\n'.format(trial_name) + 
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
            '# skip_header = 23\n' +
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
        
        for data in data_generator:
            out.write(
                '{0}     {1}     {2}     {3}     {4}     {5}     {6}     {7}     {8}     {9}\n'.format(*data)
            )

        out.close()
        
        # End

# run the monte carlo
run_sim(100, 'RUN1', output='data/images/RUN1.txt', star_min=4, verbose_factor=100)
