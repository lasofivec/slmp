#! /usr/bin/python
from scipy.sparse.linalg import spsolve, splu
from scipy.interpolate import interp2d
from scipy.io import mmread, mmwrite
import igakit.nurbs as nurbs
from scipy import interpolate
#from geometries_xml import patch
from math import pi
from post_evaluation import *
import interpol as inter
from globals_variables import *
from charac_feet import *
from distribution_functions import initialize_distribution




# *****************************************
# Defining the distribution test function

func_init = initialize_distribution(which_f)
print "=> distribution function initialized"
print "____________________________________"
print ""
# *****************************************


#============================================================#
#                                                            #
#                METHODE SEMI-LAGRANGIENNE                   #
#                                                            #
#============================================================#


advec = np.zeros((2, npatchs, NPTS1, NPTS2))

#---------------------------
#    advection definition
#---------------------------
if (which_advec == 0) :
    #    for a scalar advection :
    advec[0] = advec_dir1 #advection in 1st direction for all patches
    advec[1] = advec_dir2 #advection in 2nd direction for all patches
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



#------
zn   = np.zeros((npatchs, NPTS1, NPTS2))
znp1 = np.zeros((npatchs, NPTS1, NPTS2))
zn   = np.copy(z).reshape((npatchs, NPTS1, NPTS2))


# Plotting the initial state ***************************
plot_nrb_dens(X_mat, Y_mat, zn, func_formula, \
               show=True, save=True, tstep = -1)
# ******************************************************

# Computing the characteristics' origin
char_eta1, char_eta2, where_char = get_pat_char(geo, eta1, eta2, advec, dt)

print np.shape(char_eta1)
import numpy
XXmat = numpy.zeros_like(char_eta1)
YYmat = numpy.zeros_like(char_eta2)
for npat in range(npatchs):
    for ii in range(numpy.shape(char_eta1)[1]):
        # if numpy.shape(numpy.shape(char_eta1))[0] == 2 :
        for jj in range(numpy.shape(char_eta1)[2]):
            ee1 = char_eta1[npat][ii,jj]
            ee2 = char_eta2[npat][ii,jj]
            nnp = where_char[npat][ii,jj]
            D = geo[nnp].evaluate_deriv(ee1, ee2, nderiv=0)
            XXmat[npat,ii,jj] = D[0,0,0,0]
            YYmat[npat,ii,jj] = D[0,0,0,1]
    # else :
    #     ee1 = char_eta1[ii]
    #     ee2 = char_eta2[ii]
    #     nnp = where_char[ii]
    #     D = geo[nnp].evaluate_deriv(ee1, ee2, nderiv=0)
    #     XXmat[ii] = D[0,0,0,0]
    #     YYmat[ii] = D[0,0,0,1]

if DEBUG_MODE:
    npat = 0
    print "................", npat, "................"
    print "................ X ..............."
    print X_mat[npat]#-0.005-XXmat
    print "................ Y ..............."
    print Y_mat[npat]#-0.005-YYmat
    print "................ E1 char ..............."
    print char_eta1[npat]
    print "................ E2 char ..............."
    print char_eta2[npat]
    print "................ where ............."
    print where_char[npat]
    # print "................ Jac ..............."
    # print jac[npat]
    print "................ X adv (=F1(E1,E2))..............."
    print XXmat[npat]
    print "................ Y adv (=F1(E1,E2))..............."
    print YYmat[npat]


    
# Extracting the particles that stay in their own domain:
char_eta1_id = np.copy(char_eta1)
char_eta2_id = np.copy(char_eta2)
tab_ind_out = []
for npat in list_patchs:
    ind_out_pat = np.where(where_char[npat] != npat)
    tab_ind_out.append(ind_out_pat)
    char_eta1_id[npat][ind_out_pat] = 0.0
    char_eta2_id[npat][ind_out_pat] = 0.0

# # Using jacobian values:
# for npat in list_patchs :
#     jinv_11 =  jac[npat,3,:,:]
#     jinv_12 = -jac[npat,1,:,:]
#     jinv_21 = -jac[npat,2,:,:]
#     jinv_22 =  jac[npat,0,:,:]
#     jacob = jinv_11 * jinv_22 - jinv_21 * jinv_12
#     jacob[np.where(jacob == 0.)] = epsilon
#     char_eta1_id[npat,:,:] = eta1_mat - dt * (advec[0,0]*jinv_11 + advec[1,0]*jinv_12)/jacob
#     char_eta2_id[npat,:,:] = eta2_mat - dt * (advec[0,0]*jinv_21 + advec[1,0]*jinv_22)/jacob




# arrays for storing the errors :
list_err_inf = []
list_err_l2  = []
list_minval  = []
list_mass    = []

# ......................................................
from selalib_interpol import mp_int

for tstep in range(1,nstep+1) :

    # Computing the limit conditions :
    eta1_slopes, eta2_slopes = inter.compute_slopes(zn, list_patchs, jac)

    for npat in list_patchs :
        # Interpolation on points that stay IN the domain :
        znp1[npat] = mp_int.interpolate_2d(npat, zn[npat],
                                           char_eta1_id[npat],
                                           char_eta2_id[npat],
                                           eta1_slopes[npat],
                                           eta2_slopes[npat])
        if np.size(np.where(np.abs(znp1) > 10**5)) > 0:
            print "Gigantic value !!", np.where(np.abs(znp1) > 10**5)
            # import sys
            # sys.exit("Error !!")

        pat_out_char = tab_ind_out[npat]

        # Interpolation on points that are OUTSIDE the domain (has to be done point by point):
        for ind_pt_out in range(np.size(pat_out_char[0])):
            ind_eta1 = pat_out_char[0][ind_pt_out]
            ind_eta2 = pat_out_char[1][ind_pt_out]
            npat_char = where_char[npat][ind_eta1, ind_eta2]
            znp1[npat, ind_eta1, ind_eta2] = mp_int.f2py_interpolate_value(npat_char,
                                                                    char_eta1[npat, ind_eta1, ind_eta2],
                                                                    char_eta2[npat, ind_eta1, ind_eta2],
                                                                    zn[npat_char],
                                                                    eta1_slopes[npat_char],
                                                                    eta2_slopes[npat_char])
    zn = np.copy(znp1)
    zn[np.where(abs(zn) < 10**-10)] = 0.


    # -----------------------------------------------
    # Printing of results and time-relative error
    #------------------------------------------------
    if ((tstep == 1)or(tstep%viewstep == 0)or(tstep == nstep-1)) :
        # Plotting the density ***************************
        title  = 'Computed solution of the advection equation at $t =$ ' + str(tstep)+'\n\nwith '+func_formula
        listZ = plot_nrb_dens(X_mat, Y_mat, zn, title, show = False, save = True, tstep = tstep)
        # Computing error ***********
        # Computing the characteristic feet in PHYSICAL space at tnp1 :
        Xnp1, Ynp1 = get_phy_char(X_mat, Y_mat, npatchs, advec, which_advec, tstep, dt)
        list_errs = []
        comp_err_time(geo, X_mat, Y_mat, Xnp1, Ynp1, func_init, listZ, list_errs, plot = True, save=True, tval = tstep)
        list_err_inf.append(list_errs[0])
        list_err_l2.append(list_errs[1])
        list_minval.append(list_errs[2])
        list_mass.append(list_errs[3])
        # Printing some results
        print('= =========== TSTEP = ', tstep, "/", nstep, ' =========== =')
        maxerr = np.max(list_errs[0])
        npat_maxerr = list_errs[0].index(np.max(list_errs[0]))
        print(" --> For npat =", npat_maxerr, " maximum err l_inf =", maxerr)


# -------------------------------------------
# Computing exact solution and plotting error :
Xnp1, Ynp1 = get_phy_char(X_mat, Y_mat, npatchs, advec, which_advec, tstep, dt)
comp_err_time(geo, \
              X_mat, Y_mat, \
              Xnp1, Ynp1,
              func_init,
              listZ, list_errs,
              show = True, plot = False,
              save = True, tval = tstep)
title  = 'Computed solution of the advection equation at $t =$ ' \
         + str(tstep*dt)+'\n\nwith '+func_formula
plot_nrb_dens(X_mat, Y_mat, zn, title, show = True, save = False)

plot_errors([list_err_inf, list_err_l2, list_minval, list_mass])

##TO SAVE RESULTS :
# maxerr = np.max(final_errs[0])*10.
# file = open('resultats.txt', 'a')
# string = str(NPTS1) + " " + str(NPTS2) + " " + str(NPTS1*NPTS2) + " " + str(dt) + " " + str(tmax) + " " + str(maxerr) + "\n"
# file.write(string)
# file.close()
