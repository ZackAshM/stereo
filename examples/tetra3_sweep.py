"""
Parameter sweep to obtain solution

Goal: see which param iteration has most solves. If tie, which has lowest mean solve time.

- sweep through params for each image
- write [iteration number, fov, pme, pcs, mr, mt, solve_time] to some file
- Separately, count how many of each iteration number -> which is max? -> use those params
"""

import stereo
from stereo.sim import *
from stereo.hardware import ASI183MM_Sigma135mm as cam
import stereo.flightpath as flightpath
import tetra3
import numpy as np
import random
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from tqdm import tqdm
import os.path as op


def solution_sweep(N: int,
                   filename: str,
                   ra: float = None,
                   dec: float = None) -> None:
    """
    Sweep through FoV, pattern_match_error, pattern_checking_stars,
    match_radius, and match_threshold. See tetra3 documentation for
    descriptions of these parameters.

    Parameters
    ----------
    N : int
        The desired number of successfully solved images to achieve.
    filename : str
        The file name of the data log.
    ra, dec : float
        If given, solves for exactly one image centered at ra, dec.
    """
    if ra:
        N = 1

    # initialize parameters
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
        'pattern_checking_stars':None, 
        'match_radius':None, 
        'match_threshold':None
    }
    anita = flightpath.load_flight(4)

    solved_trials = 0

    data_dir = op.join(op.dirname(op.dirname(op.abspath(__file__))), *("data", "images"))
    output = "{0}/{1}.txt".format(data_dir, filename)
    out = open(output, 'w')
    
    while solved_trials < N:
    
        stars = 0
        # while stars < 10: # high sweep
        while 4 > stars or stars > 10: # low sweep
            # set observation
            ind = random.randrange(0, 37279)
            lat = anita.latitude[ind]
            lon = anita.longitude[ind]
            alt = anita.altitude[ind]
            time = Time(anita.realTime[ind], format="unix")
            obs = Observation(lat=lat, lon=lon, alt=alt, time=time)
            altaz = obs.altaz_frame

            # set center ra,dec
            if ra:
                center_radec = SkyCoord(ra=ra, dec=dec, unit=u.deg)
            else:
                cen_az = random.randrange(0, 360) # degrees
                center_radec = SkyCoord(alt=cen_alt, az=cen_az, unit=u.deg, frame=altaz).transform_to('icrs')

            image = generate_starfield(center_radec, cam, obs, mag_limit, psf_sigma)

            stars = len(image.image_stars)
            if stars < 4 and N == 1:
                print("Less than 4 stars in field.")
                return

        # set parameter ranges
        fovs = np.array(['db_fov6ps15cs20pme1', 'db_fov7ps15cs20pme005', 'db_fov10ps15cs20pme005'])
        pattern_match_error = np.array([0.001, 0.005, 0.01, 0.05])
        pattern_checking_stars = np.arange(4, 9)
        match_radius = np.array([0.001, 0.01])
        match_threshold = np.array([1e-10, 1e-9, 1e-8, 1e-7])

        params = [(fov, pme,pcs,mr,mt)
                  for fov in fovs
                  for pme in pattern_match_error
                  for pcs in pattern_checking_stars 
                  for mr in match_radius 
                  for mt in match_threshold]
        
        solved = False
        _fov = ''
        iteration = 0

        print("# Solutions:", solved_trials, end='\n')
            
        for fov, pme, pcs, mr, mt in tqdm(params):
            if not _fov == fov:
                db = False
            if not db:
                _fov = fov
                t3 = tetra3.Tetra3()
                t3.load_database(fov)
                db = True
            t3.database_properties['pattern_max_error'] = pme
            solve_dict['pattern_checking_stars'] = pcs
            solve_dict['match_radius'] = mr
            solve_dict['match_threshold'] = mt

            iteration += 1
            try:
                result = run_tetra3(image, t3=t3, return_result=True, **centroid_dict, **solve_dict)
            except StopIteration:
                continue

            if result['RA']:
                if np.abs(result['RA'] - center_radec.ra.value) < 6:
                    if not solved:
                        solved_trials += 1
                        solved = True
                    st = result['T_solve']
                    fovdeg = fov.split('ps')[0].split('fov')[1]
                    out.write('{0} {1} {2} {3} {4} {5} {6} {7}\n'.format(iteration, fovdeg, pme, pcs, mr, mt, st, stars))
    out.close()

# run the sweep
solution_sweep(100, filename='lowsweep2')
