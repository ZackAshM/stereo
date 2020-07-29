"""
This module provides a method to create a Camera object
with specifications for the ASI183MM CCD and Sigma 135mm F1.8 lens.
"""

import os.path as op

from stereo.sim import Camera

# define the ASI183MM + Sigma 135mm F1.8 parameters
sensor_width = 5496  # pix
sensor_height = 3672  # pix
pix = 2.4e-3  # mm/pix
focal_length = 135  # mm
fnumber = 1.8
avg_noise = 1.6
temp = 20  # C
exp_time = 0.3  # s
max_ct = 15000

# create the Camera object for ASI183MM + Sigma 135mm F1.8 setup
ASI183MM_Sigma135mm = Camera(
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

# get the directory containg the camera data
data_directory = op.join(
    op.dirname(op.dirname(op.dirname(op.abspath(__file__)))), "data"
)

# load QE and dark current
QEfile = op.join(data_directory, "ASISigma_QE.csv")
ASI183MM_Sigma135mm.loadQE(QEfile, unpack=True, delimiter=",")

dcfile = op.join(data_directory, "ASISigma_current.csv")
ASI183MM_Sigma135mm.loadDark(dcfile, unpack=True, delimiter=",")
