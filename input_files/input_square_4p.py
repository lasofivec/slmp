from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = SQUARE_4p_domain
# Discretization:
NPTS1 = 50
NPTS2 = 50

# *********************************************
# Definition of simulation
# *********************************************
tmax = 5.
dt   = 0.025
viewstep = 1

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION"
dist_center_x = 0.15
dist_center_y = 0.

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST
advec_dir1 = -0.1
advec_dir2 = 0.

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3
