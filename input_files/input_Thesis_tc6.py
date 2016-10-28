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
tmax = 1.
dt   = 0.0005
viewstep = 40

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION"
dist_center_x = 0.65
dist_center_y = 0.

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CIRC

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3
