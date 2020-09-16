"""
This module provides a function to run the tetra3 star tracking
algorithm over given generated starfield images.
"""

import os.path as op
from pathlib import Path
from typing import Any, Dict, Union

import timeout_decorator
from astropy import units as u
from astropy.coordinates import SkyCoord

import tetra3
from stereo.sim import Image


def run_tetra3(
    image: Image,
    t3: tetra3.Tetra3 = None,
    database: str = "db_fov6ps15cs20pme01",
    return_result: bool = False,
    timeout: float = 5,
    **kwargs: Any,
) -> Union[None, SkyCoord, Dict]:
    """
    Perform the tetra3 star tracking algorithm on a given image to obtain the
    estimated center sky coordinates.

    Parameters
    ----------
    image : Image
        The Image object for which the algorithm analyzes.
    t3 : tetra3.Tetra3
        If given, the Tetra3 object to use for the solver.
    database : str
        If t3 not given, the name of the database to load.
    return_result : bool
        If True, returns the direct result of tetra3.solve_from_image. Otherwise,
        returns a SkyCoord with the solution ra and dec (or None if not solved).
    timeout : float
        The time in seconds to raise a 'StopIteration' error.
    **kwargs
        Keyword arguments passed to tetra3.solve_from_image and
        tetra3.get_centroids_from_image.

    Returns
    -------
    result : Dict
        If return_result is True, the resulting dict returned from
        tetra3.solve_from_image.
    center : SkyCoord or None
        The estimated ra,dec SkyCoord of the center of the image.
        Returns None if algorithm fails to estimate a center.
    """

    # define the function using the timeout decorator
    @timeout_decorator.timeout(timeout, timeout_exception=StopIteration)
    def _run_tetra3(
        image: Image,
        t3: tetra3.Tetra3 = None,
        database: str = "db_fov6ps15cs20pme01",
        return_result: bool = False,
        **kwargs: Any,
    ) -> Union[None, SkyCoord, Dict]:

        # create Tetra3 object if not given and load a database
        if not t3:
            t3 = tetra3.Tetra3()

            data_directory = op.join(
                op.dirname(op.dirname(op.dirname(op.abspath(__file__)))),
                *("data", "tetra3"),
            )

            # load or make appropriate database
            # NOTE: To generate a database, tetra3 requires the Yale Bright Star
            # Catalog in its base directory.
            # Direct download: http://tdc-www.harvard.edu/catalogs/BSC5
            while True:
                try:
                    database_path = Path(op.join(data_directory, database))
                    t3.load_database(database_path)
                except FileNotFoundError:
                    print(
                        "Warning: No appropriate database found.",
                        "Using default database.",
                    )
                    database = "db_fov6ps15cs20pme01"
                    continue
                break

        # run the star tracking algorithm
        result = t3.solve_from_image(image.image, **kwargs)

        # return the tetra3 result
        if return_result:
            return result

        # or return the ra, dec
        if result["RA"] is None or result["Dec"] is None:
            center = None
        else:
            center = SkyCoord(ra=result["RA"], dec=result["Dec"], unit=u.deg)
        return center

    # return the result of the decorated function
    return _run_tetra3(image, t3, database, return_result, **kwargs)
