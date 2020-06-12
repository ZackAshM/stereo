from typing import Tuple

import numpy as np


def get_fov(
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
