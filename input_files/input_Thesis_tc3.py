from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = SQUARE_4p_domain
# Discretization:
NPTS1 = 100
NPTS2 = NPTS1

# *********************************************
# Definition of simulation
# *********************************************
tmax = 10.
dt   = 0.05
viewstep = 10

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION"
dist_center_x = -0.35
dist_center_y = -0.1

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST
advec_dir1 = 0.1
advec_dir2 = 0.05

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3
