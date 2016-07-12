import numpy as np
import math
from reading_geometry import get_geometry
from numpy import sin, cos, pi, sqrt

import sys
sys.path.append('./input_files/')
from input_options import *


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
which_deriv = DERIV_FDU2
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

# Getting the geometry:
NPTS    = [NPTS1, NPTS2]
geo     = get_geometry(domain)
npatchs = geo.npatchs

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

name_deriv = ""
if (which_deriv == 0) :
    name_deriv = "Appx derivative method : Finite difference centered order 2"
elif (which_deriv == 1) :
    name_deriv = "Appx derivative method : Finite difference uncentered order 2"
elif (which_deriv == 2) :
    name_deriv = "Appx derivative method : Lagrangian uncentered order 4"
elif (which_deriv == 3) :
    name_deriv = "Appx derivative method : Lagrangian centered order 4"
elif (which_deriv == 4) :
    name_deriv = "Appx derivative method : Numpy gradient function order 2"
name_deriv = name_deriv + " (" + str(which_deriv) + ")"


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
if which_f == "COS_SIN_DISTRIBUTION" :
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
     f.write('which_deriv = ' + str(which_deriv)+'\n')
     f.write('NPTS1 = '+ str(NPTS1)+'\n')
     f.write('NPTS2 = '+ str(NPTS2)+'\n')
     f.write('tmax = '+ str(tmax)+'\n')
     f.write('dt = '+ str(dt)+'\n')
     f.write('nstep = '+ str(nstep)+'\n')
     f.write('viewstep = '+ str(viewstep)+'\n')
     f.write('dist_center_x = '+ str(dist_center_x)+'\n')
     f.write('dist_center_y = '+ str(dist_center_y)+'\n')



#TODO : division of 2 integers => float ? from __future__ import division

# ***************************************
# Definition of personal mesh grid function
# ***************************************
def my_meshgrid(eta1, eta2) :
    et1, et2 = np.meshgrid(eta1, eta2)
    return et1.transpose(), et2.transpose()

epsilon = 0.01

def jacobian_function(npat, e1, e2):
    d1F1 = np.zeros((NPTS1, NPTS2))
    d2F1 = np.zeros((NPTS1, NPTS2))
    d1F2 = np.zeros((NPTS1, NPTS2))
    d2F2 = np.zeros((NPTS1, NPTS2))
    if npat == 0:
        d1F1 = 0.5*cos(-0.5*pi*e2-0.5*pi)
        d2F1 = 0.5*pi*(0.5*e1+0.5)*sin(-0.5*pi*e2-0.5*pi)
        d1F2 = 0.5*sin(-0.5*pi*e2-0.5*pi)
        d2F2 = -0.5*pi*(0.5*e1+0.5)*cos(-0.5*pi*e2-0.5*pi)
    elif npat == 1:
        d1F1 = 0.5*cos(-0.5*pi*e2)
        d2F1 = 0.5*pi*(0.5*e1+0.5)*sin(-0.5*pi*e2)
        d1F2 = 0.5*sin(-0.5*pi*e2)
        d2F2 = -0.5*pi*(0.5*e1+0.5)*cos(-0.5*pi*e2)
    elif npat == 2:
        d1F1 = -0.5*pi*(0.5*e2+0.5)*sin(0.5*pi*e1)
        d2F1 = 0.5*cos(0.5*pi*e1)
        d1F2 = 0.5*pi*(0.5*e2+0.5)*cos(0.5*pi*e1)
        d2F2 = 0.5*sin(0.5*pi*e1)
    elif npat == 3:
        d1F1 = -0.5*pi*(0.5*e2+0.5)*sin(0.5*pi*e1+0.5*pi)
        d2F1 = 0.5*cos(0.5*pi*e1+0.5*pi)
        d1F2 = 0.5*pi*(0.5*e2+0.5)*cos(0.5*pi*e1+0.5*pi)
        d2F2 = 0.5*sin(0.5*pi*e1+0.5*pi)
    elif npat == 4:
        sqrte1 = sqrt(1. - 0.5*(2.*e1-1.)**2)
        sqrte2 = sqrt(1. - 0.5*(2.*e2-1.)**2)
        e1e2 = (2.*e1-1.)*(2.*e2-1.)
        d1F1 =  0.5/sqrt(2.)*(e1e2/sqrte1 - 2.*sqrte2)
        d2F1 = -0.5/sqrt(2.)*(2.*sqrte1 - e1e2/sqrte2)
        d1F2 = -0.5/sqrt(2.)*(e1e2/sqrte1 + 2.*sqrte2)
        d2F2 =  0.5/sqrt(2.)*(2.*sqrte1 + e1e2/sqrte2)
    return [d1F1, d2F1, d1F2, d2F2]
    
