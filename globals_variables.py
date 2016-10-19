import numpy as np
import math
from reading_geometry import get_geometry
from numpy import sin, cos, pi, sqrt

import sys
sys.path.append('./input_files/')
from input_options import *


# ***************************************
# Definition of personal mesh grid function
# ***************************************
def my_meshgrid(eta1, eta2) :
    et1, et2 = np.meshgrid(eta1, eta2)
    return et1.transpose(), et2.transpose()


# Default values (see their meaning in the input_options file)
domain = SQUARE_2p_domain
NPTS1 = 48
NPTS2 = 48
tmax  = 1.2
dt    = 0.05
viewstep = 10
which_f = DIST_GAUSS
dist_center_x = 0.5
dist_center_y = 0.5
which_advec = ADVEC_CNST
which_interp = INTER_SLL3
DEBUG_MODE = False
advec_dir1 = 0.1
advec_dir2 = 0.1
# *******************************


# Reading input File ............................
from sys import argv
if len(argv) > 1:
    scriptname, filename = argv
    filename = filename[:-3]
    filename = filename.split('/')[-1]
    print "Reading the inputfile...", filename
    exec('from '+filename+' import *')
#................................................


# number of time steps
nstep = int(math.floor(tmax/dt))


name_advec = ""
if (which_advec == 0) :
    name_advec = "Constant coefficient advection"
elif (which_advec == 1) :
    name_advec = "Centered circular advection"
elif (which_advec == 2) :
    name_advec = "Uncentered circular advection"
    centera = -0.1
    centerb =  0.0
name_advec = name_advec + " (" + str(which_advec) + ")"

name_deriv = "Method to compute Hermite slopes : Cubic spline approximation"

#-------------------
# data printing
#-------------------
print ""
print " **", name_advec, " **"
print " **", name_deriv, " **"
print " NPTS1, NPTS2 =", NPTS1, NPTS2
print " dt = ", dt

#--------------------
# plotting data
#--------------------
if which_f == "CONSTANT_DISTRIBUTION" :
    PLOT_VAL_MAX = 1.
    PLOT_VAL_MIN = -1.
    comap = "RdBu" #"Spectral"
    func_formula='$f(x,y) = 0.5$'
elif which_f == "X_DISTRIBUTION" :
    PLOT_VAL_MAX = 2.
    PLOT_VAL_MIN = 0.
    comap = "RdBu"
    func_formula='$f(x,y) = x$'
elif which_f == "COS_SIN_DISTRIBUTION" :
    PLOT_VAL_MAX = 1.
    PLOT_VAL_MIN = -1.
    comap = "RdBu" #"Spectral"
    func_formula='$\cos{(2\pi x)}*\sin{(2\pi y)}$'
elif which_f == "GAUSSIAN_DISTRIBUTION":
    PLOT_VAL_MAX = 1.
    PLOT_VAL_MIN = -0.2
    comap = 'plasma'
    func_formula="$\exp(-0.5(((x+1)mod(2)-1.5)^2/0.04^2 + (y-0.5)^2/0.04^2))$"
else:
    print "Please define the proper values for vmin and vmax"
    import sys
    sys.error("Error: No defined plot_val_max nor plot_val_min")


# **************************
# function to write data
# **************************
def write_globals(path, str_num) :
     f = open(path+'/global_variables'+str_num, 'w')
     f.write("Domain = " + str(domain)+'\n')
     f.write('which_f = ' + str(which_f)+'\n')
     f.write('which_advec = ' + str(which_advec)+'\n')
     f.write('NPTS1 = '+ str(NPTS1)+'\n')
     f.write('NPTS2 = '+ str(NPTS2)+'\n')
     f.write('tmax = '+ str(tmax)+'\n')
     f.write('dt = '+ str(dt)+'\n')
     f.write('nstep = '+ str(nstep)+'\n')
     f.write('viewstep = '+ str(viewstep)+'\n')
     f.write('dist_center_x = '+ str(dist_center_x)+'\n')
     f.write('dist_center_y = '+ str(dist_center_y)+'\n')


epsilon = 0.01
epsilon2 = 10**-12
small_epsilon = 10**-14
