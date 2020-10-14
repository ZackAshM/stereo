# tetra3 files

This directory contains tetra3 databases at two different fields of view (6 and 7 degrees). The parameter 'pme' or 'pattern match error' can be changed when running a Monte Carlo sim (see examples/run_sim.py). To use them, they must be located at the top of your local tetra3 repository.

To allow generation of custom databases, the tetra3 repository requires the file 'BSC5'.

To run trail fit centroiding, it is necessary for now to replace the local copy of 'tetra3.py' in the tetra3 repository with the 'tetra3.py' file located in this directory.