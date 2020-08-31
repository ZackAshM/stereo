"""                                                                                                                                     
This module provides a function to run the tetra3 star tracking                                                                         
algorithm over given generated starfield images.                                                                                        
"""

import os.path as op
from pathlib import Path
from typing import Any, Optional

from astropy import units as u
from astropy.coordinates import SkyCoord

import tetra3
from stereo.sim import Image

import timeout_decorator


@timeout_decorator.timeout(60*5, timeout_exception=StopIteration)
def run_tetra3(image: Image,
               database: str = 'db_fov10ps15cs20pme01',
               t3: tetra3.Tetra3 = None,
               **kwargs: Any) -> Optional[SkyCoord]:
    """                                                                                                                                 
    Perform the tetra3 star tracking algorithm on a given image to obtain the                                                           
    estimated center sky coordinates.                                                                                                   
                                                                                                                                        
    Parameters                                                                                                                          
    ----------                                                                                                                          
    image : Image                                                                                                                       
        The Image object for which the algorithm analyzes.                                                                              
    **kwargs                                                                                                                            
        Keyword arguments passed to tetra3.solve_from_image and                                                                         
        tetra3.get_centroids_from_image. 

    Returns
    -------
    center : SkyCoord or None
        The estimated ra,dec SkyCoord of the center of the image.
        Returns None if algorithm fails to estimate a center.
    """

    # max_fov = image.cam.fov[2]  # diagonal
    # mag_limit = image.mag_limit

    if not t3:
        t3 = tetra3.Tetra3()

        data_directory = op.join(
            op.dirname(op.dirname(op.dirname(op.abspath(__file__)))),
            *("data", "tetra3")
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
                    "Using default database.")
                database = 'db_fov10ps15cs20pme01'
                continue
            break

    # run the star tracking algorithm
    result = t3.solve_from_image(image.image, **kwargs)

    # return the ra, dec
    if result["RA"] is None or result["Dec"] is None:
        # print("Tetra3 failed to obtain a result.")
        center = None
    else:
        center = SkyCoord(ra=result["RA"], dec=result["Dec"], unit=u.deg)

    return center
