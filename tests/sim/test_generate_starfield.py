import stereo
from stereo.sim.generate_starfield import generate_starfield

ra = 10.68
dec = 41.26
radius = 0.1
obs_time = "2016-12-26 22:38:52"
obs_lat, obs_lon, obs_alt = -85.7, 141.7, 39277
mag_limit = None

table = generate_starfield(ra, dec, radius, obs_lat, obs_lon, obs_alt, obs_time)

print(table)
