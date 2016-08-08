#! /usr/bin/python
from geometry import z, X_mat, Y_mat, jac, eta1, eta2, func_init, geo, npatchs, list_patchs
from scipy.sparse.linalg import spsolve, splu
from scipy.interpolate import interp2d
from scipy.io import mmread, mmwrite
import igakit.nurbs as nurbs
from scipy import interpolate
from math import pi
from post_evaluation import *
import interpol as inter
from globals_variables import *
from charac_feet import *
import math



#---------------------------
#    advection definition
#---------------------------
advec = np.zeros((2, npatchs, NPTS1, NPTS2))
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


#============================================================#
#                                                            #
#                METHODE SEMI-LAGRANGIENNE                   #
#                                                            #
#============================================================#


# Computing the characteristics' origin
char_eta1, char_eta2, where_char = get_pat_char(geo, eta1, eta2, advec, dt)
    
# Extracting the particles that stay in their own domain:
char_eta1_id = np.copy(char_eta1)
char_eta2_id = np.copy(char_eta2)
tab_ind_out = []
for npat in list_patchs:
    ind_out_pat = np.where(where_char[npat] != npat)
    tab_ind_out.append(ind_out_pat)
    char_eta1_id[npat][ind_out_pat] = 0.0
    char_eta2_id[npat][ind_out_pat] = 0.0


# arrays for storing the errors :
list_err_inf = []
list_err_l2  = []
list_minval  = []
list_maxval  = []
list_mass    = []

print "For pat 4 face 1 =", inter.get_face(zn[4], 1)
print "For pat 4 face 2 =", inter.get_face(zn[4], 2)


# ......................................................
from selalib_interpol import mp_int

for tstep in range(1,nstep+1) :

    # Computing the limit conditions :
    eta1_slopes, eta2_slopes = inter.compute_slopes(zn, list_patchs, jac)
    STOP
    for npat in list_patchs :
        znp1[npat] = 0.
        # Interpolation on points that stay IN the domain :
        znp1[npat] = mp_int.interpolate_2d(zn[npat],
                                           char_eta1_id[npat],
                                           char_eta2_id[npat],
                                           eta1_slopes[npat],
                                           eta2_slopes[npat])

        if np.size(np.where(np.abs(znp1) > 10**5)) > 0:
            print "Gigantic value !!", np.where(np.abs(znp1) > 10**5)
            import sys
            sys.exit("Error !!")

        if np.isnan(znp1[npat]).any() :
            print "NaN value found !!!! in Pat", npat
            print znp1[npat]
            www = np.where(np.isnan(znp1[npat]) == True)
            # print "where =", www
            #print "znp1[where] =", znp1[npat][www]
            print "zn     =", zn[npat]
            print "char 1 =", char_eta1_id[npat]
            print "char 2 =", char_eta2_id[npat]
            print "slope1 =", eta1_slopes[npat]
            print "slope2 =", eta2_slopes[npat]
            import sys
            sys.exit("Error !!")
        

        # pat_out_char = tab_ind_out[npat]

        # # Interpolation on points that are OUTSIDE the domain (has to be done point by point):
        # for ind_pt_out in range(np.size(pat_out_char[0])):
        #     ind_eta1 = pat_out_char[0][ind_pt_out]
        #     ind_eta2 = pat_out_char[1][ind_pt_out]
        #     npat_char = where_char[npat][ind_eta1, ind_eta2]
        #     znp1[npat, ind_eta1, ind_eta2] = mp_int.interpolate_value(
        #         char_eta1[npat, ind_eta1, ind_eta2],
        #         char_eta2[npat, ind_eta1, ind_eta2],
        #         zn[npat_char],
        #         eta1_slopes[npat_char],
        #         eta2_slopes[npat_char])
    zn = np.copy(znp1)
    zn[np.where(abs(zn) < 10**-10)] = 0.

    print "tstep =", tstep
    print "For pat 4 face 1 =", inter.get_face(zn[4], 1)[-3:]
    print "For pat 4 face 2 =", inter.get_face(zn[4], 2)[:3]
    print "For pat 3 face 0 =", inter.get_face(zn[3], 0)[:3]
    print "For pat 2 face 0 =", inter.get_face(zn[2], 0)[-3:]

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
        list_maxval.append(list_errs[3])
        list_mass.append(list_errs[4])
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

plot_errors([list_err_inf, list_err_l2, list_minval, list_maxval, list_mass])

##TO SAVE RESULTS :
# maxerr = np.max(final_errs[0])*10.
# file = open('resultats.txt', 'a')
# string = str(NPTS1) + " " + str(NPTS2) + " " + str(NPTS1*NPTS2) + " " + str(dt) + " " + str(tmax) + " " + str(maxerr) + "\n"
# file.write(string)
# file.close()
