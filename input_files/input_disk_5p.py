from input_options import *

DEBUG_MODE = False

# *********************************************
# Definition of the geometry
# *********************************************
domain = DISK_5p_domain
# Discretization:
NPTS1 = 40
NPTS2 = 40

# *********************************************
# Definition of simulation
# *********************************************
tmax = 15
dt   = 0.1
viewstep = 10

# *********************************************
# Definition of initial distribution
# *********************************************
which_f = "GAUSSIAN_DISTRIBUTION" 
dist_center_x = -0.65
dist_center_y = 0.2

# *********************************************
# Definition of advection coefficient
# *********************************************
which_advec = ADVEC_CNST
advec_dir1 = 0.05
advec_dir2 = -0.05

#**********************************************
# Definition of interpolation method used
# *********************************************
which_interp = INTER_SLL3 
