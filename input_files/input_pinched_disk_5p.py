from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = PINCHEDDISK_5p_domain
# Discretization:
NPTS1 = 60
NPTS2 = 60

# *********************************************
# Definition of simulation
# *********************************************
tmax = 10
dt   = 0.1
viewstep = 10

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION" 
dist_center_x = 0.2#-0.65
dist_center_y = 0.#4

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST
advec_dir1 = 0.05#-0.
advec_dir2 = 0.#-0.025

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3 
