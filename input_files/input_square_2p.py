from input_options import *

# *********************************************
# Definition of the geometry
# *********************************************
domain = SQUARE_2p_domain
# Discretization:
NPTS1 = 16
NPTS2 = 16

# *********************************************
# Definition of simulation
# *********************************************
tmax = 9.
dt   = 0.05
viewstep = 10

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION" 
dist_center_x = 0.5
dist_center_y = 0.5

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST
advec_dir1 = 0.1
advec_dir2 = 0.

#**********************************************
# Definition of derivative approximation
# *********************************************
which_deriv  = DERIV_LGC4
which_interp = INTER_SLL3 
