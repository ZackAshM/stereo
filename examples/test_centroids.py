import stereo
from stereo.sim import *
from stereo.hardware import ASI183MM_Sigma135mm as cam

from tetra3 import Tetra3, get_centroids_from_image as centroids

import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u

# define camera (ASI..Sigma can be used as cam,
# but this is to change parameters if desired)
# sensor_width = 5496 #pix
# sensor_height = 3672 #pix
# pix = 2.4e-3 #mm/pix
# focal_length = 135 #mm
# fnumber = 1.8
# avg_noise = 1.6
# temp = 20
# exp_time = 1
# max_ct = 15000
# cam1 = Camera(sensor_width,
#              sensor_height,
#              pix,
#              focal_length,
#              fnumber,
#              avg_noise,
#              temp,
#              exp_time,
#              max_ct
#             )
# cam1.loadQE('data/ASISigma_QE.csv', QE_peak=0.84, unpack=True, delimiter=',')
# cam1.loadDark('data/ASISigma_current.csv', unpack=True, delimiter=',')


# get the observation
obs = load_anita_observation(4, 0)

# define center ra,dec coords
cen_alt = 45 #degrees
cen_az = 0 # degrees
altaz = obs.altaz_frame
# center_radec = SkyCoord(alt=cen_alt, az=cen_az, unit=u.deg, frame=altaz).transform_to('icrs')
center_radec = SkyCoord(ra=84.0500, dec=-1.2019, unit=u.deg) # Orion
# center_radec = SkyCoord(ra=296.75827806346416, dec=11.314043605282833, unit=u.deg) # tetra3

mag_limit = 7
psf_sigma = 3

# generate the image
image = generate_starfield(center_radec, cam, obs, mag_limit, psf_sigma)

# get_centroids_from_image kwargs
centroid_dict = {'sigma':3, 'image_th':None, 'filtsize':10, 'binary_open':True, 'bg_sub_mode':'global_median',
                 'sigma_mode':'global_median_abs',
                 'centroid_window':None, 'max_area':None, 'min_area':5, 'max_sum':None, 'min_sum':100, 
                 'max_axis_ratio':None}

# get the centroids
star_pos = centroids(image.image, **centroid_dict)

# output the results
print(star_pos)
image.plot(centroids=star_pos)
