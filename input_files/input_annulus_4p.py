from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = ANNULUS_4p_domain
# Discretization:
NPTS1 = 51
NPTS2 = 51

# *********************************************
# Definition of simulation
# *********************************************
tmax = 8.5
dt   = 0.005
viewstep = 20

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION" 
dist_center_x = -0.8
dist_center_y = -0.15

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST

#**********************************************
# Definition of derivative approximation
# *********************************************
which_deriv  = 0#DERIV_LGC4
which_interp = INTER_SLL3 
