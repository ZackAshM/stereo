from typing import Tuple

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.table import Table
from astropy.time import Time

from scipy.ndimage import gaussian_filter


def fov(
    sensor_width: float, sensor_height: float, focal_length: float
) -> Tuple[float, float, float]:
    """
    Calculates the horizontal, vertical, and diagonal
    field of view (FoV) of a camera given its sensor
    width and height, and the focal length of the
    camera lens.

    Parameters
    ----------
    sensor_width : float
        The width of the camera sensor (mm).
    sensor_height : float
        The height of the camera sensor (mm).
    focal_length : float
        The focal length of the camera lens (mm).
    Returns
    -------
    xfov, yfov, dfov : float
        The horizontal (x), vertical (y), and
        diagonal (d) field of view in degrees.
    """

    # set the values
    w = sensor_width
    h = sensor_height
    f = focal_length
    d = np.sqrt(w ** 2 + h ** 2)

    # calculate the FoV in radians
    xfov_rad = 2 * np.arctan2(w / 2, f)
    yfov_rad = 2 * np.arctan2(h / 2, f)
    dfov_rad = 2 * np.arctan2(d / 2, f)

    # convert them to degrees
    xfov = np.rad2deg(xfov_rad)
    yfov = np.rad2deg(yfov_rad)
    dfov = np.rad2deg(dfov_rad)

    # and return them
    return xfov, yfov, dfov


def bright_conesearch(
    center_coord: SkyCoord = None, radius: float = None, mag_limit: float = None
) -> Table:
    """
    Return an astropy table of star coordinates and magnitudes
    from the Yale Bright Star Catalog. Specify a center
    sky coordinate and a radius to return the bright stars
    within a cone of the given radius from the center.
    Specify a magnitude limit to return stars whose magnitude
    are below the given limit.

    Parameters
    ----------
    center_coord : SkyCoord
        The sky coordinates of the center in a conesearch.
    radius : float
        The radius from the center_coord in a conesearch (deg).
    mag_limit : float
        The magnitude cutoff limit.

    Returns
    -------
    catalog : Table
        The catalog of bright stars. Includes only stars satisfying
        the conesearch or magnitude cutoff if they are used.
        Otherwise returns the full sky.
    """

    # get the bright star catalog
    BSC5 = "yale_bright_star_catalog5.fits.gz"
    catalog = Table.read(BSC5)

    # change column names for convenience
    catalog.rename_columns(("RAJ2000", "DEJ2000", "Vmag"), ("ra", "dec", "Mag"))

    # cutoff stars above mag_limit if given
    if mag_limit is not None:
        catalog = catalog[catalog["Mag"] <= mag_limit]

    # cutoff stars beyond radius of center if given
    if radius is not None:
        # raise error if center coord not given
        if center_coord is None:
            raise ValueError(
                "Center coordinates, center_coords=SkyCoord(), "
                "must be provided if radius is given."
            )

        # get the star coords
        cat_ra = catalog["ra"]
        cat_dec = catalog["dec"]
        starcoords = SkyCoord(ra=cat_ra, dec=cat_dec, unit=u.deg)
        catalog["skycoord"] = starcoords

        # get the separation of stars from the center
        separations = catalog["skycoord"].separation(center_coord)
        catalog["separation"] = separations

        # get only stars within the given radius
        catalog = catalog[catalog["separation"] <= radius * u.deg]

    # return the ra, dec, mag of resulting catalog
    return catalog["ra", "dec", "Mag"]


def radec2altaz(
    radec_table: Table,
    obs_lat: float,
    obs_lon: float,
    obs_alt: float,
    obs_time: str = None,
) -> Table:
    """
    Convert the (ra, dec) coordinates in a given astropy table
    of stars with columns 'ra', 'dec', and 'Mag' to (alt, az)
    using the observation location and time.

    Parameters
    ----------
    radec_table : Table
        The astropy table of stars with column coordinates
        'ra' and 'dec' in degrees, and magnitude column 'Mag'.
    obs_lat : float
        The latitude on Earth of the observation (deg).
    obs_lon : float
        The longitude on Earth of the observation (deg).
    obs_alt : float
        The altitude on Earth of the observation (m).
    obs_time : str
        The observation time (YYYY-MM-DDTHH:MM:SS).

    Returns
    -------
    altaz_table : astropy.table.Table
        The converted table with column coordinates 'alt'
        and 'az' in degrees, and magnitude column 'Mag'.
    altaz : AltAz
        The alt, az frame used for the conversion.
    """

    # make the table for alt, az coords
    altaz_table = Table(radec_table, names=["alt", "az", "Mag"], copy=True)

    # get the observation location and time
    loc = EarthLocation(lat=obs_lat * u.deg, lon=obs_lon * u.deg, height=obs_alt * u.m)
    time = Time(obs_time)

    # define the AltAz frame
    altaz = AltAz(location=loc, obstime=time)

    # get skycoords for ra,dec transformed to alt,az
    rd2aa = SkyCoord(
        ra=radec_table["ra"], dec=radec_table["dec"], obstime=obs_time, unit=u.deg
    ).transform_to(altaz)

    # set table values to the transformed coords
    altaz_table["alt"] = rd2aa.alt.value * rd2aa.alt.unit
    altaz_table["az"] = rd2aa.az.value * rd2aa.az.unit

    # return result
    return altaz_table, altaz


def map_stars(
    altaz_stars: Table,
    altaz_center: SkyCoord,
    xfov: float,
    yfov: float,
    width: int,
    height: int,
) -> np.ndarray:
    """
    Map the (alt, az) positions of stars onto pixels of a frame
    with dimensions (height, width) centered at altaz_center.

    Parameters
    ----------
    altaz_stars : Table
        Astropy table of stars with at least columns 'alt' and 'az'
        in degrees.
    altaz_center : SkyCoord
        The sky coordinates in (alt, az) of the center of the
        framewhich to map the stars onto.
    xfov : float
        The horizontal field of view in degrees.
    yfov : float
        The vertical field of view in degrees.
    width : int
        The pixel width of the frame which to map the stars onto.
    height : int
        The pixel height of the frame which to map the stars onto.

    Returns
    -------
    mapped_stars : ndarray
        An ndarray with dimensions (height, width) with pixels
        equal to 1 if a star is located there.
    """

    # get the center coordinate values
    center_alt = altaz_center.alt.value
    center_az = altaz_center.az.value

    # define the frame boundaries x=az, y=alt
    alt_lo = center_alt - (yfov / 2)
    alt_hi = center_alt + (yfov / 2)
    az_lo = center_az - (xfov / 2)
    az_hi = center_az + (xfov / 2)

    # get the alt, az arrays of the stars
    alt1 = altaz_stars["alt"]
    az1 = altaz_stars["az"]

    # remove any stars outside the boundaries
    alt_bound = np.logical_and(alt1 >= alt_lo, alt1 <= alt_hi)
    az_bound = np.logical_and(az1 >= az_lo, az1 <= az_hi)
    bound = np.logical_and(alt_bound, az_bound)
    stars = altaz_stars[bound]

    # get the alt, az of the remaining stars
    alt = stars["alt"]
    az = stars["az"]

    # create an empty frame
    frame = np.zeros((height, width))

    # get the span of alt, az positions
    azspan = az_hi - az_lo
    altspan = alt_hi - alt_lo

    # get the mapping scale
    azscale = (az - az_lo) / azspan
    altscale = (alt - alt_lo) / altspan

    # map the positions to the frame indices
    azind = np.round_(azscale * (width - 1)).astype(int)
    altind = np.round_(altscale * (height - 1)).astype(int)

    # set pixel values of alt, az in the frame to 1
    frame[altind, azind] = 1

    # and return the frame
    return frame


def generate_starfield(
    center_radec: SkyCoord,
    cam_params: dict,
    obs_params: dict,
    mag_limit: float = 7,
    psf_sigma: float = 1,
) -> np.ndarray:
    """
    Generate a star field centered around center_coord
    (ra,dec) using camera parameters and observation parameters.
    A cutoff magnitude can be specified if desired. Stars follow
    a Gaussian Point Spread Function with a specified standard
    deviation.

    Parameters
    ----------
    center_radec : SkyCoord
        The (ra, dec) sky coordinates of the center of the frame.
    cam_params : Dict[str, float]
        Camera parameters 'sensor_width' (mm), 'sensor_height' (mm),
        and 'focal_length' (mm).
    obs_params : Dict[str, float]
        Observation parameters 'obs_lat' (deg), 'obs_lon' (deg),
        'obs_alt' (deg), and 'obs_time' (str).
    mag_limit : float
        All stars used will be of magnitude less than this.
    psf_sigma : float
        The size/standard deviation of the Guassian Point Spread Function.

    Returns
    -------
    starfield : ndarray
        The array of pixel values making up the star field.
    """

    # calculate the field of view
    xfov, yfov, dfov = fov(**cam_params)

    # generate a table of bright stars within fov
    radius = dfov / 2
    stars_radec = bright_conesearch(center_radec, radius, mag_limit)

    # convert the star coordinates to alt, az
    stars_altaz, altaz = radec2altaz(stars_radec, **obs_params)

    # get the center alt, az
    center_altaz = center_radec.transform_to(altaz)

    # map the remaining stars to a frame
    width = cam_params["sensor_width"]
    height = cam_params["sensor_height"]
    starframe = map_stars(stars_altaz, center_altaz, xfov, yfov, width, height)

    # apply a Gaussian PSF
    starfield = gaussian_filter(starframe, sigma=psf_sigma)

    # and return the starfield
    return starfield
