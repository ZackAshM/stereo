"""
This module provides the Image class.
"""

from typing import Dict, Tuple

import attr
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from numpy import logical_and as AND
from numpy import logical_not as NOT
from scipy.integrate import simps
from scipy.ndimage import gaussian_filter

from stereo.sim.camera import Camera
from stereo.sim.observation import Observation


@attr.s
class Image(object):
    """
    An image class containing at least a table of stars
    and a sky center, and capable of image generation of stars
    using Camera and Observation objects.

    Attributes
    ----------
    stars : astropy.table.Table
        The table of stars 'ra', 'dec', and 'Mag'. Not necessarily the
        stars shown within the image.
    center : astropy.coordinates.SkyCoord
        The sky coordinates of the center of the image.
    observation : Observation
        The Observation object.
    cam : Camera
        The Camera object.
    psf_sigma : float
        The standard deviation of the Gaussian PSF in pixels. Default = 2.

    Properties
    ----------
    size : Tuple[int, int]
        The size of the image frame in pixels.
    altaz_stars : astropy.table.Table
        The table of stars 'alt', 'az', 'Mag'.
    image_bounds : Dict[str, Tuple[float, float]]
        The alt, az bounds of the image.
    image_stars : astropy.table.Table
        The table of stars 'alt', 'az', 'Mag' shown in the image.
    image : numpy.ndarray
        The ndarray of a pixel grid containing pixel values
        corresponding to the location and magnitudes of the stars.

    Methods
    -------
    mag2pix(table: Table) : astropy.table.Table
        Returns the ra, dec stars table or a given table
        with magnitude values converted to total pixel counts.
    plot()
        Minimally plots the image data.
    """

    stars: Table = attr.ib()
    center: SkyCoord = attr.ib()
    observation: Observation = attr.ib(default=None)
    cam: Camera = attr.ib(default=None)
    psf_sigma: float = attr.ib(default=2)

    @property
    def size(self) -> Tuple[int, int]:
        """Returns the pixel dimensions of the CCD."""
        try:
            return self._size  # type: ignore
        except AttributeError:
            if self.cam is None:
                raise ValueError("Camera object required to determine image size")
            size = (self.cam.sensor_pixheight, self.cam.sensor_pixwidth)
            self._size = size
            return size

    @property
    def altaz_stars(self) -> Table:
        """
        Convert the (ra, dec) coordinates in a given astropy table
        of stars with columns 'ra', 'dec', and 'Mag' to (alt, az)
        using the observation location and time.

        Returns
        -------
        altaz_table : astropy.table.Table
            The converted table with column coordinates 'alt'
            and 'az' in degrees, and magnitude column 'Mag'.

        Raises
        ------
        ValueError
            If Observation object attribute is not set.
        """

        try:
            return self._altaz_stars  # type: ignore
        except AttributeError:
            if self.observation is None:
                raise ValueError("Observation object required for alt/az conversion")

            # make the table for alt, az coords
            altaz_table = Table(self.stars, names=["alt", "az", "Mag"], copy=True)

            # get skycoords for ra,dec transformed to alt,az
            obs_time = self.observation.time
            rd2aa = SkyCoord(
                ra=self.stars["ra"], dec=self.stars["dec"], obstime=obs_time, unit=u.deg
            ).transform_to(self.observation.altaz_frame)

            # set table values to the transformed coords
            altaz_table["alt"] = rd2aa.alt.value * rd2aa.alt.unit
            altaz_table["az"] = rd2aa.az.value * rd2aa.az.unit

            # return result
            self._altaz_stars = altaz_table
            return altaz_table

    @property
    def image_bounds(self) -> Dict[str, Tuple[float, float]]:
        """
        Returns the alt, az bounds of the image.

        Returns
        -------
        bounds : Dict[str, Tuple[float, float]]
            The dictionary containing the 'alt' and 'az' bounds as (lower, upper)
            tuples.

        Raises
        ------
        ValueError
            If Camera or Observation object attributions are not set.
        """
        try:
            return self._image_bounds  # type: ignore
        except AttributeError:
            if None in [self.cam, self.observation]:
                raise ValueError(
                    "Camera and Observation objects required to"
                    + " determine image bounds"
                )

            # use alt, az stars
            altaz_center = self.center.transform_to(self.observation.altaz_frame)
            center_alt = altaz_center.alt.value
            center_az = altaz_center.az.value

            # define the frame boundaries x=az, y=alt
            xfov, yfov, _ = self.cam.fov
            alt_lo = center_alt - (yfov / 2)
            alt_hi = center_alt + (yfov / 2)
            az_lo = center_az - (xfov / 2)
            az_hi = center_az + (xfov / 2)

            bounds = {"alt": (alt_lo, alt_hi), "az": (az_lo, az_hi)}

            # return the boundaries
            self._image_bounds = bounds
            return bounds

    @property
    def image_stars(self) -> Table:
        """Returns the 'alt', 'az', and 'Mag' table of stars shown in the image."""
        try:
            return self._image_stars  # type: ignore
        except AttributeError:

            # get alt, az and bounds
            star_table = Table(self.altaz_stars, copy=True)

            alt = star_table["alt"]
            az = star_table["az"]

            alt_lo, alt_hi = self.image_bounds["alt"]
            az_lo, az_hi = self.image_bounds["az"]

            # circles means modulo 360
            alt_lo = alt_lo % 360
            alt_hi = alt_hi % 360
            az_lo = az_lo % 360
            az_hi = az_hi % 360
            alt = alt % 360
            az = az % 360

            # determine boolean bounds
            if alt_lo < alt_hi:
                alt_bound = AND(alt_lo <= alt, alt <= alt_hi)
            else:
                alt_bound = NOT(AND(alt_hi <= alt, alt <= alt_lo))

            if az_lo < az_hi:
                az_bound = AND(az_lo <= az, az <= az_hi)
            else:
                az_bound = NOT(AND(az_hi <= az, az <= az_lo))

            # remove stars outside the bounds
            bound = np.logical_and(alt_bound, az_bound)
            stars = star_table[bound]

            # and return them
            self._image_stars = stars
            return stars

    @property
    def image(self) -> np.ndarray:
        """
        Returns the ndarray containing pixel values corresponding to
        positions and magnitudes in self.stars with a Gaussian PSF.

        Returns
        -------
        image_data : ndarray
            An ndarray with dimensions of CCD with pixels
            equal to the total pixel flux where a star is located.
        """

        try:
            return self._image  # type: ignore
        except AttributeError:

            image_stars = self.image_stars
            stars = self.mag2pix(image_stars)

            alt = stars["alt"]
            az = stars["az"]
            ct = stars["Total Counts"]

            height, width = self.size
            image_data = np.zeros((height, width))

            # get the span of alt, az positions
            alt_lo, alt_hi = self.image_bounds["alt"]
            az_lo, az_hi = self.image_bounds["az"]
            azspan = az_hi - az_lo
            altspan = alt_hi - alt_lo

            # get the mapping scale
            azscale = (az - az_lo) / azspan
            altscale = (alt - alt_lo) / altspan

            # map the positions to the frame indices
            azind = np.round_(azscale * (width - 1)).astype(int)
            altind = np.round_(altscale * (height - 1)).astype(int)

            # set pixel values at star pixel positions
            for i in zip(altind, azind, ct):
                image_data[i[0], i[1]] = i[2]

            # apply Gaussian PSF
            image_data = gaussian_filter(image_data, sigma=self.psf_sigma)

            # add noise and dark current
            image_data += np.random.poisson(self.cam.avg_noise, size=self.size)
            noise_level = np.interp(self.cam.temp, *self.cam.dark_current)
            image_data += np.random.lognormal(
                noise_level * self.cam.exp_time, size=self.size
            )

            image_data = np.clip(image_data, 0, self.cam.max_ct)

            # and return the image
            self._image = image_data
            return image_data

    def mag2pix(self, table: Table = None) -> Table:
        """
        Returns the stars table with magnitudes converted to total pixel counts.

        Parameters
        ----------
        table : astropy.table.Table
            The star table to use. If None, the ra,dec stars table is used.

        Returns
        -------
        stars : Table
            The table containing pixel values corresponding to the magnitudes.
        """

        # get the desired star table
        if table is not None:
            stars = Table(table, copy=True)
        else:
            stars = Table(self.stars, copy=True)

        # quantify the magnitudes
        mag = np.array(stars["Mag"])
        mag_val = mag * u.STmag

        # convert magnitudes to flux
        flux = mag_val.to(
            u.photon / u.s / u.cm ** 2 / u.nm, u.spectral_density(550 * u.nm)
        )

        # get the total QE response
        QE = self.cam.QE
        if QE is not None:
            response = simps(QE[1], QE[0]) * u.nm / u.photon
        else:
            response = 600 * u.nm / u.photon  # 100% from 400nm to 1000nm

        aper_area = np.pi * (self.cam.radius * 0.1) ** 2

        # calculate the pixel counts
        pixel_cts = flux * self.cam.exp_time * aper_area * response

        # replace the magnitudes with pixel values
        stars["Mag"] = pixel_cts
        stars.rename_column("Mag", "Total Counts")

        # and return the table
        return stars

    def plot(self) -> None:
        """Minimally plot the image data."""
        fig, ax = plt.subplots(figsize=[20, 10])
        plt.gray()
        plt.title(
            r"ra{0:.2f}$\degree$dec{1:.2f}$\degree$, FoV: {2:.1f}$\degree$".format(
                self.center.ra.value, self.center.dec.value, self.cam.fov[2]
            )
        )
        ax.imshow(self.image, norm=LogNorm(vmin=1, vmax=self.cam.max_ct))
        plt.show()
