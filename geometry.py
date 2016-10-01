import numpy as np
from globals_variables import dist_center_x
from globals_variables import dist_center_y
from globals_variables import NPTS1, NPTS2
from globals_variables import domain
from globals_variables import my_meshgrid
from globals_variables import which_f
from reading_geometry import get_geometry
from reading_geometry import get_patches
from distribution_functions import initialize_distribution


# Getting the geometry:
NPTS    = [NPTS1, NPTS2]
geo     = get_geometry(domain)
npatchs = geo.npatchs
list_patchs = get_patches(geo)

# *****************************************
# Defining the distribution test function

func_init = initialize_distribution(which_f,
                                    center_x=dist_center_x,
                                    center_y=dist_center_y)
print "=> distribution function initialized"
print "____________________________________"
print ""
# *****************************************

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
