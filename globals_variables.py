import numpy as np
import math
from reading_geometry import get_geometry
from reading_geometry import get_patches
from distribution_functions import initialize_distribution
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


# *****************************************
# Defining the distribution test function

func_init = initialize_distribution(which_f, dist_center_x=dist_center_x, dist_center_y=dist_center_y)
print "=> distribution function initialized"
print "____________________________________"
print ""
# *****************************************

# ***************************************
# Definition of personal mesh grid function
def my_meshgrid(eta1, eta2) :
    et1, et2 = np.meshgrid(eta1, eta2)
    return et1.transpose(), et2.transpose()
# ***************************************


# number of time steps
nstep = int(math.floor(tmax/dt))

# Getting the geometry:
NPTS    = [NPTS1, NPTS2]
geo     = get_geometry(domain)
npatchs = geo.npatchs
list_patchs = get_patches(geo)

# Computing the physical coo & jacobian of the transformation:
#-------------------------------------------
# Defining knots and value of field on knots
#-------------------------------------------
eta1  = np.linspace(0., 1., NPTS1)
eta2  = np.linspace(0., 1., NPTS2)
z     = np.zeros((npatchs, NPTS1*NPTS2))
X_mat = np.zeros((npatchs, NPTS1, NPTS2))
Y_mat = np.zeros((npatchs, NPTS1, NPTS2))
jac   = np.zeros((npatchs, 4, NPTS1, NPTS2))

eta1_mat, eta2_mat = my_meshgrid(eta1,eta2)
# *********************************************

for npat in list_patchs :
    # Calculating the corresponding values
    # of knots on the physical space :
    D = geo[npat].evaluate_deriv(eta1,eta2, nderiv=1)
    X_mat[npat] = D[0,:,:,0]
    Y_mat[npat] = D[0,:,:,1]
    # Calculation the density on these points :
    X       = X_mat[npat].reshape((NPTS1*NPTS2))
    Y       = Y_mat[npat].reshape((NPTS1*NPTS2))
    z[npat] = func_init(X, Y)
    # python has an UF for very small values of x,y at exp(x,y) :
    z[np.where(abs(z) < 10**-9)] = 0.
    # Computing jacobian values
    jac[npat,0,:,:] = D[1, :,:, 0]
    jac[npat,1,:,:] = D[2, :,:, 0]
    jac[npat,2,:,:] = D[1, :,:, 1]
    jac[npat,3,:,:] = D[2, :,:, 1]


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


epsilon = 0.01
epsilon2 = 10**-14

