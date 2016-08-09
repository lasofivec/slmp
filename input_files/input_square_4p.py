from input_options import *

# *********************************************
# Definition of the geometry
# *********************************************
domain = SQUARE_4p_domain
# Discretization:
NPTS1 = 100
NPTS2 = 100

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
dist_center_x = -0.5
dist_center_y = -0.5

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3 
