from input_options import *

# *********************************************
# Definition of the geometry
# *********************************************
domain = ANNULUS_4p_domain
# Discretization:
NPTS1 = 101
NPTS2 = 101

# *********************************************
# Definition of simulation
# *********************************************
tmax = 7.5
dt   = 0.05
viewstep = 10

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
which_deriv  = DERIV_LGC4
which_interp = INTER_SLL3 
