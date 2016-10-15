from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = PINCHEDDISK_5p_domain
# Discretization:
NPTS1 = 100
NPTS2 = 100

# *********************************************
# Definition of simulation
# *********************************************
tmax = 9
dt   = 0.05
viewstep = 50

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION"
dist_center_x = -0.65
dist_center_y = 0.

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST
advec_dir1 = 0.15
advec_dir2 = 0.

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3
