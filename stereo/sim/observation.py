"""
This module provides the Observation class.
"""

import attr
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation
from astropy.time import Time


@attr.s
class Observation(object):
    """
    The parameters describing an observation location and time.

    Attributes
    ----------
    lat : float
        The latitude in degrees.
    lon : float
        The longitude in degrees.
    alt : float
        The altitude in meters.
    time : str
        The time of observation in any astropy acceptable time formats.

    Properties
    ----------
    altaz_frame : astropy.coordinates.AltAz
        The frame needed to transform ra, dec to alt, az.
    """

    lat: float = attr.ib()
    lon: float = attr.ib()
    alt: float = attr.ib()
    time: Time = attr.ib(default=Time("2000-01-01 00:00:00", format="iso"))

    @property
    def altaz_frame(self) -> AltAz:
        """Returns the frame for transforming ra,dec to alt,az."""
        try:
            return self._altaz_frame  # type: ignore
        except AttributeError:
            loc = EarthLocation(
                lat=self.lat * u.deg, lon=self.lon * u.deg, height=self.alt * u.m
            )
            time = Time(self.time)
            altaz = AltAz(location=loc, obstime=time)
            self._altaz_frame = altaz
            return altaz
