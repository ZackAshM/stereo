"""
This module provides methods to load the flightpath
data of various ANITA flights.
"""
import collections
import os.path as op
from typing import Any

import uproot

# the directory where we store our data
data_directory = op.join(op.dirname(op.dirname(__file__)), "data")


def load_flight(flight: int = 4) -> Any:
    """
    Load the flightpath for a given version of ANITA.

    The returned array has *at least* the following fields:

        realTime:  the time of each entry in unix time.
        altitude:  the payload altitude in m.
        latitude:  the payload latitude in degrees.
        longitude: the payload longitude in degrees.
        heading:   the payload heading in degrees.

    For ANITA3 and ANITA4, it also contains

        pitch: the payload pitch in degrees.
        roll:  the payload roll in degrees.

    Parameters
    ----------
    flight: int
        The flight number to load.

    Returns
    -------
    flightpath: uproot.tree.Arrays
        A namedtuple-like class containing numpy arrays for each quantity.
    """

    # check for a valid version
    if flight not in [1, 3, 4, 5]:
        raise ValueError(f"We currently only support ANITA{1, 3, 4, 5} (got: {flight})")

    # if we are simulating, use ANITA-4's flight path
    if flight == 5:
        flight = 4

    # construct the filename for this flight
    filename = op.join(data_directory, *("flightpaths", f"anita{flight}.root"))

    # open the ROOT file
    f = uproot.open(filename)

    # and load the ttree and return it to the user
    return f["adu5PatTree"].arrays(outputtype=collections.namedtuple)
