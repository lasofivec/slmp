from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = DISK_5p_domain
# Discretization:
NPTS1 = 100
NPTS2 = 100

# *********************************************
# Definition of simulation
# *********************************************
tmax = 1.5
dt   = 0.005
viewstep = 10

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION"
dist_center_x = 0.62
dist_center_y = 0.4

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CIRC

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3
