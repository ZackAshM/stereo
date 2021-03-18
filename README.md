# stereo

[![Actions Status](https://github.com/rprechelt/stereo/workflows/Pytest/badge.svg)](https://github.com/rprechelt/stereo/actions)
![GitHub](https://img.shields.io/github/license/rprechelt/stereo?logoColor=brightgreen)
![Python](https://img.shields.io/badge/python-3.6%20%7C%203.7%20%7C%203.8-blue)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A daytime star tracker for the Payload for Ultra-high Energy Observations (PUEO)

A package to perform Monte Carlo simulations to investigate the performance of an open-source daytime star tracking algorithm, Tetra3, on the determination of the in-flight attitude of the PUEO payload. The MC generates low-magnitude star fields using Yale's Bright Star Catalogue, with variable signal-to-noise ratio and flight rotation speed, by randomly choosing a sky center, earth location, etc. based off of the ANITA (Antarctive Impulsive Transient Antenna) flight data. Results show the performance of attitude determination for varying SNR and rotation speeds, as well as different centroiding algorithms (averaging vs trail fitting).
