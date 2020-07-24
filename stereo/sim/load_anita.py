from astropy.time import Time

import stereo.flightpath as flightpath
from stereo.sim import Observation


def load_anita_observation(flight: int = 4, ind: int = 0) -> Observation:
    """
    Load event latitude, longitude, altitude, and time from
    ANITA flight data into an Observation object.

    Parameters
    ----------
    flight : {1, 2, 3, 4}
        The flight number. Default = 4.
    ind : int
        The index to determine the event used. Default = 0.

    Returns
    -------
    obs: Observation
        The Observation object holding the event data: latitude (deg),
        longitude (deg), altitude (m), and time.
    """

    # arrays: heading, latitude, longitude, altitude, realTime, pitch, roll
    data = flightpath.load_flight(flight)

    lat = data[1][ind]
    lon = data[2][ind]
    alt = data[3][ind]
    time = Time(data[4][ind], format="unix")

    obs = Observation(lat=lat, lon=lon, alt=alt, time=time)

    return obs
