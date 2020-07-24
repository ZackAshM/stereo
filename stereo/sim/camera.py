"""
This module provides the classes Camera, Observation, and Image.
"""

from typing import Any, Tuple

import attr
import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.table import Table
from astropy.time import Time
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from scipy.integrate import simps
from scipy.ndimage import gaussian_filter


@attr.s
class Camera(object):
    """
    Camera parameters used for star trackers. Includes information
    on both the CCD and the lens.

    Attributes
    ----------
    sensor_pixwidth : int
        The CCD width in pixels.
    sensor_pixheight : int
        The CCD height in pixels.
    pix_size : float
        The size of CCD pixels in mm/pix.
    focal_length : float
        The lens focal length in mm.
    fnumber : float
        The lens f-number.
    avg_noise : float
        The average noise of the CCD. Default = 0.
    temp : float
        The temperature of the CCD in celsius. Default = 20.
    exp_time : float
        The exposure time in seconds. Default = 0.3.
    max_ct : int
        The maximum well depth of the CCD in counts. Default = 15000.
    QE : ndarray
        The response vs wavelength (nm) of the CCD
        quantum efficiency.
    dark_current : ndarray
        The dark current (e/s) vs temperature (C) of the CCD.

    Properties
    ----------
    radius : float
        The first order radius of the lens aperture.
    fov : Tuple[float, float, float]
        The horizontal, vertical, and diagonal fields of view.

    Methods
    -------
    loadQE(path, **kwargs)
        Load QE data into the QE attribute using np.genfromtxt.
    loadDark(path, **kwargs)
        Load dark current data into the dark_current attribute
        using np.genfromtxt.
    """

    sensor_pixwidth: int = attr.ib(default=None)
    sensor_pixheight: int = attr.ib(default=None)
    pix_size: float = attr.ib(default=None)
    focal_length: float = attr.ib(default=None)
    fnumber: float = attr.ib(default=None)
    avg_noise: float = attr.ib(default=0.0)
    temp: float = attr.ib(default=20.0)
    exp_time: float = attr.ib(default=1.0)
    max_ct: int = attr.ib(default=15000)
    QE: np.ndarray = attr.ib(default=None)
    dark_current: np.ndarray = attr.ib(default=None)

    @property
    def radius(self) -> float:
        """Returns the aperture radius."""
        try:
            return self._radius  # type: ignore
        except AttributeError:
            if None in [self.focal_length, self.fnumber]:
                raise TypeError("non-None values required for focal_length and fnumber")
            radius = self.focal_length / self.fnumber / 2
            self._radius = radius
            return radius

    @property
    def fov(self) -> Tuple[float, float, float]:
        """
        Returns the horizontal, vertical, and diagonal
        field of view (FoV) in degrees.
        """
        try:
            return self._fov  # type: ignore
        except AttributeError:
            # set the values
            w = self.sensor_pixwidth
            h = self.sensor_pixheight
            f = self.focal_length
            pix = self.pix_size
            if None in [w, h, f, pix]:
                raise TypeError(
                    "non-None values required for sensor_pixwidth, ",
                    "sensor_pixheight, focal_length, and pix_size",
                )
            f = f / pix
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
            self._fov = (xfov, yfov, dfov)
            return (xfov, yfov, dfov)

    def loadQE(self, path: str, **kwargs: Any) -> None:
        """
        Load data into QE.

        Parameters
        ----------
        path: str
            The data path.
        **kwargs: Any
            Pass kwargs into numpy.genfromtxt.
        """
        self.QE = np.genfromtxt(path, **kwargs)

    def loadDark(self, path: str, **kwargs: Any) -> None:
        """
        Load data into dark_current.

        Parameters
        ----------
        path: str
            The data path.
        **kwargs: Any
            Pass kwargs into numpy.genfromtxt.
        """
        self.dark_current = np.genfromtxt(path, **kwargs)


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


@attr.s
class Image(object):
    """
    An image class containing at least a table of stars
    and a sky center, and capable of image generation of stars
    using Camera and Observation objects.

    Attributes
    ----------
    stars : astropy.table.Table
        The table of stars 'ra', 'dec', and 'Mag'.
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
    image_data : numpy.ndarray
        The ndarray of a pixel grid containing pixel values
        corresponding to the location and magnitudes of the stars.

    Methods
    -------
    mag2pix(coords) : astropy.table.Table
        Returns the stars table (of either ra,dec or alt,az)
        with magnitude values converted to total pixel counts.
    plot(**kwargs)
        Plots the image data.
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

            if None in [self.cam, self.observation]:
                raise ValueError("Camera and Observation objects required for image")

            altaz_center = self.center.transform_to(self.observation.altaz_frame)
            center_alt = altaz_center.alt.value
            center_az = altaz_center.az.value

            # define the frame boundaries x=az, y=alt
            xfov, yfov, _ = self.cam.fov
            alt_lo = center_alt - (yfov / 2)
            alt_hi = center_alt + (yfov / 2)
            az_lo = center_az - (xfov / 2)
            az_hi = center_az + (xfov / 2)

            star_table = self.mag2pix(self.altaz_stars)
            alt1 = star_table["alt"]
            az1 = star_table["az"]

            # remove any stars outside the boundaries
            alt_bound = np.logical_and(alt1 >= alt_lo, alt1 <= alt_hi)
            az_bound = np.logical_and(az1 >= az_lo, az1 <= az_hi)
            bound = np.logical_and(alt_bound, az_bound)
            stars = star_table[bound]

            alt = stars["alt"]
            az = stars["az"]
            ct = stars["Total Counts"]

            height, width = self.size
            image_data = np.zeros((height, width))

            # get the span of alt, az positions
            azspan = az_hi - az_lo
            altspan = alt_hi - alt_lo

            # get the mapping scale
            azscale = (az - az_lo) / azspan
            altscale = (alt - alt_lo) / altspan

            # map the positions to the frame indices
            azind = np.round_(azscale * (width - 1)).astype(int)
            altind = np.round_(altscale * (height - 1)).astype(int)

            # set pixel values
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
