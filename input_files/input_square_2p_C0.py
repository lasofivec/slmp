from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = SQUARE_2p_C0_domain
# Discretization:
NPTS1 = 100
NPTS2 = 100

# *********************************************
# Definition of simulation
# *********************************************
tmax = 10.
dt   = 0.75
viewstep = 1

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION" 
dist_center_x = 1.5
dist_center_y = 0.65

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST
advec_dir1 = -0.1
advec_dir2 = -0.025

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3
