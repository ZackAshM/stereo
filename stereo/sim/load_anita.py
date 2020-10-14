"""
This module provides a method to load ANITA data as an Observation object.
"""

from astropy.time import Time

import stereo.flightpath as flightpath
from stereo.sim.observation import Observation


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

    maxind = data[0].size - 1

    while True:
        try:
            lat = data.latitude[ind]
            lon = data.longitude[ind]
            alt = data.altitude[ind]
            time = Time(data.realTime[ind], format="unix")
            break
        except IndexError:
            print(
                "WARNING: the requested event index, {},".format(ind),
                "is outside the range of ANITA data of size {}.".format(maxind),
                "Defaulting to event index 0.",
            )
            ind = 0

    obs = Observation(lat=lat, lon=lon, alt=alt, time=time)

    return obs
