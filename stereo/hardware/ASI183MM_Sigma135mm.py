"""
This module provides a method to create a Camera object
with specifications for the ASI183MM CCD and Sigma 135mm F1.8 lens.
"""

from stereo.sim import Camera


def ASI183MM_Sigma135mm() -> Camera:
    """
    Returns a Camera object with specifications for the
    ASI183MM CCD and Sigma 135mm F1.8 lens.
    """

    sensor_width = 5496  # pix
    sensor_height = 3672  # pix
    pix = 2.4e-3  # mm/pix
    focal_length = 135  # mm
    fnumber = 1.8
    avg_noise = 1.6
    temp = 20
    exp_time = 1
    max_ct = 15000

    camera = Camera(
        sensor_width,
        sensor_height,
        pix,
        focal_length,
        fnumber,
        avg_noise,
        temp,
        exp_time,
        max_ct,
    )

    camera.loadQE("data/ASISigma_QE.csv", unpack=True, delimiter=",")
    camera.loadDark("data/ASISigma_current.csv", unpack=True, delimiter=",")

    return camera
