# The purpose of this file is to define all the multiple choice
# variables that are in the input files.
# The program will still run without any input data,
# it will take all the default values (shown here with a * and
# in globals_variables.py file). 


DEBUG_MODE = False

# Domain:
SQUARE_2p_domain        = 0 # 2 squares domain (2 patches, C1)   [*DEFAULT*]
CIRCLE_1p_domain        = 1 # circle domain (1 patch)
DISK_5p_domain          = 2 # crown (4 patches) + internal disk  ( = 5 patches)
PINCHEDDISK_5p_domain   = 3 # pinched 5 patch disk (5 patches, no singular points)
SQUARE_4p_domain        = 4 # 2 squares domain (2 patches, C1)
ANNULUS_4p_domain       = 5 # annulus divided in 4 patches (like the disk_5p but without interior patch)
SQUARE_2p_C0_domain     = 6 # 2 square not C1 only C0

# Discretization:
NPTS1 = 24 # Number of points along x
NPTS2 = 24 # Number of points along y

# Simulation:
tmax = 1.2    # Maximum time for simulation
dt   = 0.05   # Delta_t = time step
viewstep = 10 # Frequency of file and result saving

# Initial distribution:
DIST_CONST = "CONSTANT_DISTRIBUTION" # constant distribution
DIST_GAUSS = "GAUSSIAN_DISTRIBUTION" # gaussian pulse [*DEFAULT*]
DIST_COSSN = "COS_SIN_DISTRIBUTION"  # cos(x)*sin(y)
DIST_SINX  = 2 # sinusoidal waves along x-axis
DIST_SINY  = 3 # sinusoidal waves along y-axis
DIST_SINXY = 4 # sinusoidal pavement
DIST_PETRI = 5 # Petri

# Coordinates of where the initial distribution is centered:
# (these values should change with the domain, if not sure better left undefined)
dist_center_x = 0.5 # For SQUARE_2p_domain
dist_center_y = 0.5 # For SQUARE_2p_domain

# Definition of advection coefficient:
ADVEC_CNST = 0 # constant coefficient advection
#advec_dir1 = 0.1  #advection in 1st direction for all patches
#advec_dir2 = 0.1  #advection in 2nd direction for all patches
ADVEC_CIRC = 1 # circular advection (centered at 0,0)
ADVEC_UCIR = 2 # un centered circular advection
# centera = -0.1 # For ADVEC_UCIR
# centerb = 0.0 # For ADVEC_UCIR


# Method of derivative approximation:
# Determines the formula to calculate the derivatives used
# for the BC. So far the best results have been obtained used which_deriv = 1
DERIV_FDC2 = 0 # Finite Difference method Centered of order 2
DERIV_FDU2 = 1 # Finite Difference method Uncentered of order 2
DERIV_LGU4 = 2 # Lagrangian Uncentered method of order 4
DERIV_LGC4 = 3 # Lagrangian Centered method of order 4
DERIV_NPG2 = 4 # Numpy Gradient function of order 2

# Method of interpolation:
INTER_SLL3 = 0 # selalib cubic spline method
INTER_SCI1 = 1 # python's scipy linear  spline interpolator
INTER_SCI3 = 2 # python's scipy cubic   spline interpolator
INTER_SCI5 = 3 # python's scipy quantic spline interpolator

