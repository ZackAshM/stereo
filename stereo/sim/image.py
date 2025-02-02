"""
This module provides the Image class.
"""

from typing import Any, Dict, Tuple

import attr
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
from numpy import logical_and as AND
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
        The ra, dec sky coordinates of the center of the image.
    observation : Observation
        The Observation object.
    cam : Camera
        The Camera object.
    mag_limit : float
        The dimmest magnitude cutoff. Default = max magnitude in stars['Mag'].
    psf_sigma : float
        The standard deviation of the Gaussian PSF in pixels. Default = 3.
    snr : float
        If specified, forces the image to have this signal-to-noise ratio
        defined by mean(star count) / sqrt(mean(noise count)).
    rotation : float
        The camera rotation in arcsec/sec.
    roll : float
        The payload roll in radians. If not given, it's the mean roll from
        ANITAIV flight data.

    Properties
    ----------
    size : Tuple[int, int]
        The size of the image frame in pixels.
    trail_length : int
        The length of star trail in pixels due to rotation.
    altaz_stars : astropy.table.Table
        The table of stars in 'alt', 'az', 'Mag'.
    proj_stars : astropy.table.Table
        The table of stars projected in 'x', 'y', 'Mag'.
    image_bounds : Dict[str, Tuple[float, float]]
        The projected x, y bounds of the image.
    image_stars : astropy.table.Table
        The table of stars 'x', 'y', 'Mag' within the image bounds.
    image : numpy.ndarray
        The ndarray of a pixel grid containing pixel values
        corresponding to the location and magnitudes of the stars.

    Methods
    -------
    radec2proj(radec: SkyCoord = None,
               ra: float = None,
               dec: float = None) : Tuple[float, float]
        Converts given ra, dec coordinates to image projection coordinates.
    mag2pix(table: Table) : astropy.table.Table
        Returns the ra, dec stars table or a given table
        with magnitude values converted to total pixel counts.
    plot(vmin: float = 1,
         vmax: float = None,
         centroids: ndarray = None,
         centroid_radius: float = 20,
         filename: str = None)
        Plots the image. Use vmin, vmax to set log range. If
        vmax is None, vmax is set to the camera's max pixel count.
        If centroids is given, plots circle annotations at given
        centroid locations. Saves the plot as filename if not None.
    """

    stars: Table = attr.ib()
    center: SkyCoord = attr.ib()
    observation: Observation = attr.ib(default=None)
    cam: Camera = attr.ib(default=None)
    mag_limit: float = attr.ib()
    psf_sigma: float = attr.ib(default=3)
    snr: float = attr.ib(default=None)
    rotation: float = attr.ib(default=0)
    roll: float = attr.ib(default=0.20560889)

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
    def trail_length(self) -> int:
        """Returns the trail length of stars in the image in pixels."""
        try:
            return self._trail_length  # type: ignore
        except AttributeError:
            if self.cam is None:
                raise ValueError("Camera object required to determine trail length")
            trail_length = abs(
                int(self.rotation * self.cam.exp_time / self.cam.pixscale)
            )
            trail_length += 1
            self._trail_length = trail_length
            return trail_length

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
        coordinates x, y.

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
            ra, dec = self.stars["ra"], self.stars["dec"]

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
                raise ValueError("Camera object required to determine image bounds")

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
        """Returns the 'x', 'y', and 'Mag' table of stars within the image bounds."""
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
        positions and magnitudes in self.stars assuming a Gaussian PSF.

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

            x_lo, x_hi = self.image_bounds["x"]
            y_lo, y_hi = self.image_bounds["y"]

            # map to pixels
            xind = np.interp(x, (x_lo, x_hi), (0.5, width)).astype(int)
            yind = np.interp(y, (y_lo, y_hi), (0.5, height)).astype(int)

            # set pixel values at star pixel positions
            for i in zip(yind, xind, ct):
                for shift in range(self.trail_length):
                    try:
                        # Coords generated at the end of trail
                        yshift = int(i[0] + shift * np.sin(self.roll))
                        xshift = int(i[1] + shift * np.cos(self.roll))
                        # Use below to have coords generated at center of trail
                        # yshift = int(i[0] + (shift - 0.5 * self.trail_length) *
                        #              np.sin(self.roll))
                        # xshift = int(i[1] + (shift - 0.5 * self.trail_length) *
                        #              np.cos(self.roll))

                        image_data[yshift, xshift] = i[2] / self.trail_length
                    except IndexError:
                        continue

            # apply Gaussian PSF
            image_data = gaussian_filter(image_data, sigma=self.psf_sigma)

            # generate a background
            read_noise = np.random.poisson(self.cam.avg_noise, size=self.size)
            dark_current_level = np.interp(self.cam.temp, *self.cam.dark_current)
            dark_current = np.random.lognormal(
                dark_current_level * self.cam.exp_time, size=self.size
            )
            background = read_noise + dark_current

            # force signal-to-noise ratio if desired
            signal = np.mean(image_data[image_data > 0])
            noise = np.mean(background[background > 0]) ** 0.5
            if self.snr is not None:
                image_data *= self.snr / (signal / noise)
            else:
                self.snr = signal / noise

            # add background
            image_data += background

            # rescale for nice plot behavior
            if np.max(image_data):
                ct_scale = self.cam.max_ct / np.max(image_data)
                image_data *= ct_scale

            # abide by well depth if necessary
            image_data = np.clip(image_data, 0, self.cam.max_ct)

            # correct reflections
            image_data = image_data[::-1, ::-1]

            # and return the image
            self._image = image_data
            return image_data

    def radec2proj(
        self, radec: SkyCoord = None, ra: float = None, dec: float = None
    ) -> Tuple[float, float]:
        """
        Convert ra, dec coordinates to the image projection
        coordinates in units per focal length. Coordinates
        can be given separately or as one SkyCoord object.
        """

        # get center and star ra, dec coords
        cen_ra, cen_dec = (
            self.center.ra.value * (np.pi / 180),
            self.center.dec.value * (np.pi / 180),
        )

        if radec is not None:
            ra, dec = radec.ra.value, radec.dec.value
        elif ra is None:
            raise ValueError(
                "Class Image method radec2proj was given no ra or dec coordinates"
            )

        ra = ra * (np.pi / 180)
        dec = dec * (np.pi / 180)

        # calculate projection coords
        # https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Book%3A_Celestial_Mechanics_(Tatum)/11%3A_Photographic_Astrometry/11.02%3A_Standard_Coordinates_and_Plate_Constants
        x = np.sin(ra - cen_ra) / (
            np.sin(cen_dec) * np.tan(dec) + np.cos(cen_dec) * np.cos(ra - cen_ra)
        )
        y = (np.tan(dec) - np.tan(cen_dec) * np.cos(ra - cen_ra)) / (
            np.tan(cen_dec) * np.tan(dec) + np.cos(ra - cen_ra)
        )

        # and return the new coordinates
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

        if self.cam is None:
            raise ValueError("Camera object required to determine pixel counts")

        # get the desired star table
        if table is not None:
            stars = Table(table, copy=True)
        else:
            stars = Table(self.stars, copy=True)

        # convert magnitudes to flux
        mag = np.array(stars["Mag"])  # Vmag
        wav_V = 551 * u.nm
        flux_V0 = (3640 * u.Jy).to(
            u.photon / u.nm / u.cm / u.cm / u.s, equivalencies=u.spectral_density(wav_V)
        )
        Vflux = 10.0 ** (-0.4 * mag)
        flux = flux_V0 * Vflux

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
        vmin: float = None,
        vmax: float = None,
        centroids: np.ndarray = None,
        centroid_radius: float = 20,
        filename: str = None,
        **subplots_kws: Any,
    ) -> None:
        """
        Plot the image.

        Parameters
        ----------
        vmin, vmax : float
            The min and max passed into the pyplot LogNorm scale.
            Default = image min/max.
        centroids : ndarray
            The centroids of the image. If not None, returns the image with position
            annotations.
        centroid_radius : float
            The radius of the centroid annotations.
        filename : str
            If given, saves the image with this filename.
            Extension is assumed .png.
        **subplots_kws
            Keyword arguments passed to plt.subplots
        """

        if "figsize" not in subplots_kws:
            subplots_kws["figsize"] = [20, 10]
        fig, ax = plt.subplots(**subplots_kws)
        plt.gray()

        if self.snr:
            title = (
                r"ra{0:.2f}$\degree$dec{1:.2f}$\degree$, "
                + r"FoV: {2:.1f}$\degree$, SNR: {3:.1f}, $\omega$: "
                + r'{4:0}"/s'
            ).format(
                self.center.ra.value,
                self.center.dec.value,
                self.cam.fov[2],
                self.snr,
                self.rotation,
            )
        else:
            title = (
                r"ra{0:.2f}$\degree$dec{1:.2f}$\degree$, "
                + r'FoV: {2:.1f}$\degree$, $\omega$: {3:0}"/s'
            ).format(
                self.center.ra.value,
                self.center.dec.value,
                self.cam.fov[2],
                self.rotation,
            )
        plt.title(title)
        vmin = np.min(self.image) if vmin is None else vmin
        vmax = np.max(self.image) if vmax is None else vmax
        ax.imshow(self.image, norm=LogNorm(vmin=vmin, vmax=vmax))
        ax.set_axis_off()

        # plot centroid annotations if desired
        if centroids is not None:
            x = centroids[:, 1]
            y = centroids[:, 0]
            for xx, yy in zip(x, y):
                circ = Circle((xx, yy), centroid_radius, ec="red", fill=False)
                ax.add_patch(circ)

        # save if desired
        if filename:
            plt.savefig(filename + ".png", bbox_inches="tight", pad_inches=0)

        plt.show()
