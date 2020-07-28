import stereo
from stereo.sim import *
from stereo.hardware import ASI183MM_Sigma135mm as cam

from astropy.coordinates import SkyCoord
from astropy import units as u

import cProfile, pstats, io
pr = cProfile.Profile()
pr.enable()


obs = load_anita_observation(4, 0)

cen_alt = 45 #degrees
cen_az = 0 # degrees
altaz = obs.altaz_frame
center_radec = SkyCoord(alt=cen_alt, az=cen_az, unit=u.deg, frame=altaz).transform_to('icrs')

mag_limit = 7
psf_sigma = 3

image = generate_starfield(center_radec, cam, obs, mag_limit, psf_sigma)

foo = image.image


pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
# pr.dump_stats('profile.txt')
