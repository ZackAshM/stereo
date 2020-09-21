"""
This module provides the Camera class.
"""

from typing import Any, Tuple

import attr
import numpy as np


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
    pixscale : float
        The ccd pixel scale in arcsec/pixel.

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

    @property
    def pixscale(self) -> float:
        try:
            return self._pixscale # type: ignore
        except AttributeError:
            fov = self.fov[0] * 60 * 60 # in arcsec
            pix = self.sensor_pixwidth
            pixscale = fov / pix
            self._pixscale = pixscale
            return pixscale

    def loadQE(self, path: str, QE_peak: float = 1.0, **kwargs: Any) -> None:
        """
        Load data into QE.

        Parameters
        ----------
        path: str
            The data path.
        QE_peak : float
            The peak QE to be specified if data is relative.
        **kwargs: Any
            Pass kwargs into numpy.genfromtxt.
        """
        QE = np.genfromtxt(path, **kwargs)
        QE[1] *= QE_peak
        self.QE = QE

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
