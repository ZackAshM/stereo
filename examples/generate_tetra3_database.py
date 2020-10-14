"""
Generate a tetra3 database with the listed parameters below.
"""

from tetra3 import Tetra3

t3 = Tetra3()

fov_db = 10 #deg
pme = 0.005 # error in ratio from pattern edges over longest edge
ps = 15 # max stars to use for patterns
cs = 20 # max stars to be detected in fov
mag = 10

db_name = "db_fov{0:0}ps{1:0}cs{2:0}pme".format(fov_db, ps, cs) + str(pme)[2:]

db_kwargs = {'max_fov':fov_db, 
            'save_as':db_name, 
            'pattern_stars_per_fov':ps,
            'catalog_stars_per_fov':cs, 
            'star_min_magnitude':mag,
            'pattern_max_error':pme}

t3.generate_database(**db_kwargs)
