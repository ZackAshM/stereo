# Uses ANITA data (angle data) to calculate the field of view

import stereo
from stereo.sim.get_fov import get_fov

sw, sh, fl = 40, 20, 135

print(get_fov(sw, sh, fl))
