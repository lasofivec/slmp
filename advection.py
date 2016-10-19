import numpy as np
from geometry import npatchs
from globals_variables import which_advec
from globals_variables import NPTS1, NPTS2
from globals_variables import advec_dir1, advec_dir2

#---------------------------
#    advection definition
#---------------------------
advec = np.zeros((2, npatchs, NPTS1, NPTS2))
if (which_advec == 0) :
    #    for a scalar advection :
    advec[0] = advec_dir1
    advec[1] = advec_dir2
    print " Constant advection coefficients: "
    for ipat in range(0, npatchs):
        print "     For patch", ipat, " advection vector: ", advec[:, ipat, 0, 0]
    print "___________________________________"
    print ""
if (which_advec == 1) :
    #    for a centered circular motion advection : ########
    advec = [-2.*np.pi*Y_mat, 2.*np.pi*X_mat]
if (which_advec == 2) :
    #    for an un centered circular motion advection : ########
    advec = [Y_mat - centerb, -X_mat + centera]
