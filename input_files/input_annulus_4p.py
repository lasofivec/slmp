from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = ANNULUS_4p_domain
# Discretization:
NPTS1 = 100
NPTS2 = 100

# *********************************************
# Definition of simulation
# *********************************************
tmax = 25.
dt   = 0.2
viewstep = 10

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION" 
dist_center_x = -0.65
dist_center_y = 0.4

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST
advec_dir1 = -0.
advec_dir2 = -0.025

#**********************************************
# Definition of derivative approximation
# *********************************************
which_deriv  = 0#DERIV_LGC4
which_interp = INTER_SLL3 
