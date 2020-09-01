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
from matplotlib.patches import Circle
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
    mag_limit : float
        The most dim magnitude cutoff. Default = max magnitude in stars['Mag'].
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
    plot(vmin: float = 1, vmax: float = None)
        Minimally plots the image data. Use vmin, vmax to set log range. If
        vmax is None, vmax is set to the camera's max pixel count.
    """

    stars: Table = attr.ib()
    center: SkyCoord = attr.ib()
    observation: Observation = attr.ib(default=None)
    cam: Camera = attr.ib(default=None)
    mag_limit: float = attr.ib()
    psf_sigma: float = attr.ib(default=2)

    @mag_limit.default
    def mag_limit_default(self) -> float:
        """Sets default mag_limit to dimmest magnitude given."""
        return np.amax(self.stars["Mag"])

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
    def proj_stars(self) -> Table:
        """
        Convert the (ra, dec) coordinates in a given astropy table
        of stars with columns 'ra', 'dec', and 'Mag' to image projected
        coordinates using the center ra, dec.

        Returns
        -------
        proj_table : astropy.table.Table
            The converted table with column coordinates 'y' and 'x' in 
            units per focal length, and magnitude column 'Mag'.
        """

        try:
            return self._proj_stars  # type: ignore
        except AttributeError:

            # make the table for projected coords
            proj_table = Table(self.stars, names=["y", "x", "Mag"], copy=True)
            
            # get star ra, dec coords
            ra, dec = self.stars['ra'], self.stars['dec']
            
            # get the projected coordinates
            x, y = self.radec2proj(ra=ra, dec=dec)

            # set table values to the transformed coords
            proj_table["y"] = y / u.deg
            proj_table["x"] = x / u.deg

            # return result
            self._proj_stars = proj_table
            return proj_table


    @property
    def image_bounds(self) -> Dict[str, Tuple[float, float]]:
        """
        Returns the projected x, y bounds of the image.

        Returns
        -------
        bounds : Dict[str, Tuple[float, float]]
            The dictionary containing the 'x' and 'y' bounds as (lower, upper)
            tuples.

        Raises
        ------
        ValueError
            If Camera object attribution is not set.
        """
        try:
            return self._image_bounds  # type: ignore
        except AttributeError:
            if self.cam is None:
                raise ValueError(
                    "Camera object required to determine image bounds"
                )
                
            cen_x, cen_y = self.radec2proj(radec=self.center)

            xfov, yfov, _ = self.cam.fov
            x_lo = cen_x - np.tan(xfov * (np.pi / 180) / 2)
            x_hi = cen_x + np.tan(xfov * (np.pi / 180) / 2)
            y_lo = cen_y - np.tan(yfov * (np.pi / 180) / 2)
            y_hi = cen_y + np.tan(yfov * (np.pi / 180) / 2)

            bounds = {"x": (x_lo, x_hi), "y": (y_lo, y_hi)}

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
            star_table = Table(self.proj_stars, copy=True)

            x = star_table["x"]
            y = star_table["y"]

            x_lo, x_hi = self.image_bounds["x"]
            y_lo, y_hi = self.image_bounds["y"]
            
            # determine boolean bounds
            x_bound = AND(x_lo <= x, x <= x_hi)
            y_bound = AND(y_lo <= y, y <= y_hi)

            # remove stars outside the bounds
            bound = np.logical_and(y_bound, x_bound)
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

            x = stars["x"]
            y = stars["y"]
            ct = stars["Total Counts"]

            height, width = self.size
            image_data = np.zeros((height, width))

            # get the span of alt, az positions
            x_lo, x_hi = self.image_bounds["x"]
            y_lo, y_hi = self.image_bounds["y"]
            
            # map to pixels
            xind = np.interp(x, (x_lo, x_hi), (0.5, width+0.5)).astype(int)
            yind = np.interp(y, (y_lo, y_hi), (0.5, height+0.5)).astype(int)

            # set pixel values at star pixel positions
            for i in zip(yind, xind, ct):
                image_data[i[0], i[1]] = i[2]

            # apply Gaussian PSF
            image_data = gaussian_filter(image_data, sigma=self.psf_sigma)
            if np.max(image_data):
                ct_scale = self.cam.max_ct / np.max(image_data)
                image_data *= ct_scale

            # add noise and dark current
            image_data += np.random.poisson(self.cam.avg_noise, size=self.size)
            noise_level = np.interp(self.cam.temp, *self.cam.dark_current)
            image_data += np.random.lognormal(
                noise_level * self.cam.exp_time, size=self.size
            )

            image_data = np.clip(image_data, 0, self.cam.max_ct)

            # flip horizontally
            image_data = image_data[:, ::-1]

            # and return the image
            self._image = image_data
            return image_data
        
    def radec2proj(self, radec: SkyCoord = None, ra: float = None, dec: float = None) -> Tuple[float, float]:
        """Convert ra, dec coordinates to the image projection coordinates in units per focal length."""
        
        # get center and star ra, dec coords
        cen_ra, cen_dec = self.center.ra.value*(np.pi / 180), self.center.dec.value*(np.pi / 180)
        if ra is None and dec is None:
            ra, dec = radec.ra.value, radec.dec.value
            
        ra = ra * (np.pi / 180)
        dec = dec * (np.pi / 180)
        
        # calculate projection coords:
        # https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Book%3A_Celestial_Mechanics_(Tatum)/11%3A_Photographic_Astrometry/11.02%3A_Standard_Coordinates_and_Plate_Constants
        x = np.sin(ra - cen_ra) / ( np.sin(cen_dec) * np.tan(dec) 
                                   + np.cos(cen_dec) * np.cos(ra - cen_ra) )
        y = ( np.tan(dec) - np.tan(cen_dec) * np.cos(ra - cen_ra) ) / ( 
            np.tan(cen_dec) * np.tan(dec) + np.cos(ra - cen_ra) )
        
        return (x, y)

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
        mag_val = mag * u.ABmag

        # convert magnitudes to flux
        flux = mag_val.to(
            u.photon / u.s / u.cm ** 2 / u.nm, u.spectral_density(551 * u.nm)
        )

        # get the total QE response
        QE = self.cam.QE
        if QE is not None:
            response = simps(QE[1], QE[0]) * u.nm / u.photon
        else:
            response = 600 * u.nm / u.photon  # 100% from 400nm to 1000nm

        # note: radius in mm -> cm
        aper_area = np.pi * (self.cam.radius * 0.1) ** 2

        # calculate the pixel counts
        pixel_cts = flux * self.cam.exp_time * aper_area * response

        # replace the magnitudes with pixel values
        stars["Mag"] = pixel_cts
        stars.rename_column("Mag", "Total Counts")

        # and return the table
        return stars

    def plot(
        self,
        vmin: float = 1,
        vmax: float = None,
        centroids: np.ndarray = None,
        centroid_radius: float = 20,
        save: bool = False,
        filename: str = None,
    ) -> None:
        """
        Plot the image.

        Parameters
        ----------
        vmin, vmax : float
            The min and max passed into the pyplot LogNorm scale.
        centroids : ndarray
            The centroids of the image. If not None, returns the image with position
            annotations.
        centroid_radius : float
            The radius of the centroid annotations.
        save : bool
            Save the image as a png if true.
        filename : str
            The save filename if save is True. Extension is assumed .png. Required
            if save is True.
        """

        fig, ax = plt.subplots(figsize=[20, 10])
        plt.gray()
        plt.title(
            r"ra{0:.2f}$\degree$dec{1:.2f}$\degree$, FoV: {2:.1f}$\degree$".format(
                self.center.ra.value, self.center.dec.value, self.cam.fov[2]
            )
        )
        vmax = self.cam.max_ct if vmax is None else vmax
        ax.imshow(self.image, norm=LogNorm(vmin=vmin, vmax=vmax))

        if centroids is not None:
            x = centroids[:, 1]
            y = centroids[:, 0]
            for xx, yy in zip(x, y):
                circ = Circle((xx, yy), centroid_radius, ec="red", fill=False)
                ax.add_patch(circ)

        if save:
            if filename is None:
                raise ValueError('When saving image, "filename" argument is required.')
            plt.savefig(filename + ".png")

        plt.show()
