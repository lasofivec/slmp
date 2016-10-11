from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = SQUARE_2p_domain
# Discretization:
NPTS1 = 30
NPTS2 = NPTS1

# *********************************************
# Definition of simulation
# *********************************************
tmax = 11.
dt   = 0.05
viewstep = 10

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "COS_SIN_DISTRIBUTION"#"X_DISTRIBUTION"

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST
advec_dir1 = 0.1
advec_dir2 = 0.0

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3
