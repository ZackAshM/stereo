"""
This module provides the methods for star field data generation,
including a conesearch method for bright stars.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from stereo.sim.camera import Camera
from stereo.sim.image import Image
from stereo.sim.observation import Observation


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
    BSC5 = "data/yale_bright_star_catalog5.fits.gz"
    catalog = Table.read(BSC5)

    catalog.rename_columns(("RAJ2000", "DEJ2000", "Vmag"), ("ra", "dec", "Mag"))

    # cutoff stars above mag_limit if given
    if mag_limit is not None:
        catalog = catalog[catalog["Mag"] <= mag_limit]

    # cutoff stars beyond radius of center if given
    if radius is not None:
        if center_coord is None:
            raise ValueError(
                "Center coordinates, center_coords=SkyCoord(), "
                "must be provided if radius is given."
            )

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


def generate_starfield(
    center_radec: SkyCoord,
    cam: Camera,
    obs: Observation,
    mag_limit: float = 7,
    psf_sigma: float = 3,
) -> np.ndarray:
    """
    Generate a star field centered around center_coord
    (ra,dec) using camera parameters and observation parameters.
    A cutoff magnitude can be specified if desired. Stars follow
    a Gaussian Point Spread Function with a specified standard
    deviation in pixels.

    Parameters
    ----------
    center_radec : SkyCoord
        The ra,dec sky coordinates of the center of the frame.
    cam : Camera
        Camera onject describing the camera properties.
    obs : Observation
        Observation object describing location and time of observation.
    mag_limit : float
        All stars used will be of magnitude less than this.
    psf_sigma : float
        The size/standard deviation of the Guassian Point Spread Function.

    Returns
    -------
    starfield : Image
        The Image object containing the star field data.
    """

    # calculate the field of view
    xfov, yfov, dfov = cam.fov

    # generate a table of bright stars within fov
    radius = dfov / 2
    stars_radec = bright_conesearch(center_radec, radius, mag_limit)

    # construct the image
    starfield = Image(
        stars=stars_radec,
        center=center_radec,
        observation=obs,
        cam=cam,
        psf_sigma=psf_sigma,
    )

    # and return the starfield
    return starfield
