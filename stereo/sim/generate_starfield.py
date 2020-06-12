# import matplotlib.pyplot as plt
# import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.table import Table
from astropy.time import Time

# import astroquery
# from astroquery.simbad import Simbad
from astroquery.vo_conesearch.conesearch import conesearch


def generate_starfield(
    ra: float,
    dec: float,
    radius: float,
    obs_lat: float,
    obs_lon: float,
    obs_alt: float,
    obs_time: str,
    mag_limit: float = None,
) -> Table:  # Temporary return
    """
    Generate the star field of a region centered at the
    sky coordinates (ra, dec) with given radius at a
    given time.

    Parameters
    ----------
    ra : float
        The RA sky coord of the center of the starfield (deg).
    dec : float
        The DEC sky coord of the center of the starfield (deg).
    radius : float
        The radius of the starfield region (deg).
    obs_lat : float
        The latitude on Earth of the observation (deg).
    obs_lon : float
        The longitude on Earth of the observation (deg).
    obs_alt : float
        The altitude on Earth of the observation (m).
    obs_time : str
        The observation time (YYYY-MM-DDTHH:MM:SS).
    mag_limit : float
        An upper limit on the star magnitudes, if desired.

    Returns
    -------
    starfield : ndarray
        The ndarray of pixel and pixel values of the
        starfield image generated.
    """

    # get the center skycoord
    # WARNING!!! This search may not do well for large radii.
    coords = SkyCoord(ra=ra, dec=dec, obstime=obs_time, unit=u.deg)

    # perform a conesearch centered at the coord above
    cone = conesearch(coords, radius * u.deg)

    # get all object of interest's ra, dec, and mag
    # limit the magnitude if desired
    if mag_limit is not None:
        objects_radec = cone["ra", "dec", "Mag"][cone["Mag"] < mag_limit]
    else:
        objects_radec = cone["ra", "dec", "Mag"]

    # get the observation location and time
    loc = EarthLocation(lat=obs_lat * u.deg, lon=obs_lon * u.deg, height=obs_alt * u.m)
    time = Time(obs_time)

    # define the AltAz frame
    altaz = AltAz(location=loc, obstime=time)

    # make another table for alt, az coords
    objects_altaz = Table(objects_radec, names=["alt", "az", "Mag"], copy=True)

    # get skycoords for ra,dec transformed to alt,az
    rd2aa = SkyCoord(
        ra=objects_radec["ra"], dec=objects_radec["dec"], obstime=obs_time, unit=u.deg
    ).transform_to(altaz)

    # set table values to the transformed coords
    objects_altaz["alt"] = rd2aa.alt.value * rd2aa.alt.unit
    objects_altaz["az"] = rd2aa.az.value * rd2aa.az.unit

    # remove stars not within rectangular region

    # place Gaussian at each alt, az

    # return result
    return objects_altaz  # temporary return
