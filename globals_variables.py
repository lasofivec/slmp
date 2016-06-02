import numpy as np
import math
from reading_geometry import get_geometry

# **********************
# Useful math function
# **********************
from numpy import sin, cos, pi, sqrt

# ***************************************
# Definition of the geometry
# ***************************************
SQUARE_2p_domain = 0 #--> 2 squares domain (2 patches, C1)
CIRCLE_1p_domain = 1 #--> circle domain (1 patch)
DISK_5p_domain = 2 #--> crown (4 patches) + internal disk  ( = 5 patches)
PINCHEDDISK_5p_domain = 3 #--> pinched 5 patch disk (5 patches, no singular points ?)
SQUARE_4p_domain = 4 #--> 2 squares domain (2 patches, C1)


domain = DISK_5p_domain
# Scalars for the deformation :
li_scal_x = 1.0
li_scal_y = 1.0
# Getting the geometry:
geo         = get_geometry(domain)
npatchs     = geo.npatchs
#****************************************


# *********************************************
# Definition of initial distribution
# *********************************************
# which_f = "CONSTANT_DISTRIBUTION" #--> constant distribution
which_f = "GAUSSIAN_DISTRIBUTION" #--> gaussian pulse
# which_f = "COS_SIN_DISTRIBUTION" # cos(x)*sin(y)
# which_f = 2 #--> sinusoidal waves along x-axis
# which_f = 3 #--> sinusoidal waves along y-axis
# which_f = 4 #--> sinusoidal pavement
# which_f = 5 #--> Petri
xmax = 1.0
ymax = 1.0
#****************************************

# *********************************************
# Definition of advection coefficient
# *********************************************
# which_advec = 0 --> constant coefficient advection
# which_advec = 1 --> circular advection (centered at 0,0)
# which_advec = 2 --> un centered circular advection
# which_advec = 3 --> #
which_advec = 0
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


center_adv_x = 0.25
center_adv_y = 0.25
if (domain == SQUARE_2p_domain) :
    center_adv_x = 0.5
    center_adv_y = 0.5
if (domain == SQUARE_4p_domain) :
    center_adv_x = -0.5
    center_adv_y = -0.5
if (domain == CIRCLE_1p_domain) :
    center_adv_x = -0.5
    center_adv_y = -0.25
if domain == DISK_5p_domain :
    center_adv_x = -0.65
    center_adv_y = -0.38
# *********************************************


#**********************************************
# Definition of derivative approximation
# *********************************************
# which_deriv determines the formula to calculate the derivatives used
# for the BC. So far the best results have been obtained used which_deriv = 1
# expected values 0, 1, 2 or 3
# 0,1,2 : finite diff, changes the order of apprx
# 3 : lagrangian of order 2
which_deriv = 3
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

# *** To choose interpolation function : ***
# expected values 0, 1, 2, 3
# if 0 then uses selalib cubic spline method
# if 1 then uses python's scipy linear  spline interpolator
# if 2 then uses python's scipy cubic   spline interpolator
# if 3 then uses python's scipy quantic spline interpolator
which_interp = 0

# *********************************************
# data
# *********************************************
NPTS1 = 3
NPTS2 = 3

NPTS  = [NPTS1, NPTS2]
tmax = 9.
dt = 0.05
if which_advec == 1:
    tmax = 1.0
    dt = 0.004
nstep = int(math.floor(tmax/dt))
viewstep = 10


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


#-----------------
# circle :
#-----------------
rayon = np.sqrt(center_adv_x**2 + center_adv_y**2)
anglec = np.linspace(0., 2.*np.pi, 100)
circlex =  rayon * np.cos(anglec)
circley =  rayon * np.sin(anglec)



# **************************
# function to write data
# **************************
def write_globals(path, str_num) :
     f = open(path+'/global_variables'+str_num, 'w')
     f.write("Domain = " + str(domain)+'\n')
     f.write('li_scal_x = ' + str(li_scal_x)+'\n') 
     f.write('li_scal_y = ' + str(li_scal_y)+'\n')
     f.write('which_f = ' + str(which_f)+'\n')
     f.write('xmax = ' + str(xmax)+'\n')
     f.write('ymax = ' + str(ymax)+'\n')
     f.write('which_advec = ' + str(which_advec)+'\n')
     f.write('which_deriv = ' + str(which_deriv)+'\n')
     f.write('NPTS1 = '+ str(NPTS1)+'\n')
     f.write('NPTS2 = '+ str(NPTS2)+'\n')
     f.write('tmax = '+ str(tmax)+'\n')
     f.write('dt = '+ str(dt)+'\n')
     f.write('nstep = '+ str(nstep)+'\n')
     f.write('viewstep = '+ str(viewstep)+'\n')
     f.write('center_adv_x = '+ str(center_adv_x)+'\n')
     f.write('center_adv_y = '+ str(center_adv_y)+'\n')



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
    
