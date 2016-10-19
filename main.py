#! /usr/bin/python
from geometry import z, X_mat, Y_mat, jac, eta1, eta2, npatchs, list_patchs
from scipy.sparse.linalg import spsolve, splu
from scipy.interpolate import interp2d
from scipy.io import mmread, mmwrite
import igakit.nurbs as nurbs
from scipy import interpolate
from math import pi
from post_evaluation import *
import connectivity as conn
import derivatives as deriv
from globals_variables import *
from advection import advec
from charac_feet import *
from selalib_interpol import mp_int
import math


# arrays for storing the errors :
list_err_inf = []
list_err_l2  = []
list_minval  = []
list_maxval  = []
list_mass    = []

# arrays for storing the function :
zn   = np.zeros((npatchs, NPTS1, NPTS2))
znp1 = np.zeros((npatchs, NPTS1, NPTS2))
zn   = np.copy(z).reshape((npatchs, NPTS1, NPTS2))


# Plotting the initial state ***************************
plot_nrb_dens(zn, show=True, save=True, tstep = 0)
# ******************************************************


# Computing the characteristics' origin
char_eta1, char_eta2, where_char = get_pat_char(eta1, eta2, advec, dt)


# Extracting the particles that stay in their own domain:
char_eta1_id = np.copy(char_eta1)
char_eta2_id = np.copy(char_eta2)
tab_ind_out = []
for npat in list_patchs:
    ind_out_pat = np.where(where_char[npat] != npat)
    tab_ind_out.append(ind_out_pat)
    char_eta1_id[npat][ind_out_pat] = 0.0
    char_eta2_id[npat][ind_out_pat] = 0.0


# ........................ TIME LOOP ..............................
for tstep in range(1,nstep+1) :

    # Computing the limit conditions :
    eta1_slopes, eta2_slopes = deriv.compute_slopes(zn, list_patchs, jac)

    for npat in list_patchs :
        znp1[npat] = 0.

        # Interpolation on points that stay IN the domain :
        znp1[npat] = mp_int.interpolate_2d(zn[npat],
                                           char_eta1_id[npat],
                                           char_eta2_id[npat],
                                           eta1_slopes[npat],
                                           eta2_slopes[npat])

        # Interpolation on points that are OUTSIDE the domain
        # (has to be done point by point):
        pat_out_char = tab_ind_out[npat]
        for ind_pt_out in range(np.size(pat_out_char[0])):
            ind_eta1 = pat_out_char[0][ind_pt_out]
            ind_eta2 = pat_out_char[1][ind_pt_out]
            npat_char = where_char[npat][ind_eta1, ind_eta2]
            znp1[npat, ind_eta1, ind_eta2] = mp_int.interpolate_value(
                char_eta1[npat, ind_eta1, ind_eta2],
                char_eta2[npat, ind_eta1, ind_eta2],
                zn[npat_char],
                eta1_slopes[npat_char],
                eta2_slopes[npat_char])

        # Some control test
        if np.size(np.where(np.abs(znp1[npat]) > 10**5)) > 0:
            print "Gigantic value !! in Pat ", npat
            import sys
            sys.exit("Error !!")

        if np.isnan(znp1[npat]).any() :
            print "NaN value found !!!! in Pat ", npat
            import sys
            sys.exit("Error !!")


    zn = np.copy(znp1)
    zn[np.where(abs(zn) < 10**-10)] = 0.

    # -----------------------------------------------
    # Printing of results and time-relative error
    #------------------------------------------------
    if ((tstep == 1)or(tstep%viewstep == 0)or(tstep == nstep-1)) :
        # Plotting the density:
        listZ = plot_nrb_dens(zn, tstep = tstep)
        # Computing the characteristic feet in PHYSICAL space at tnp1 :
        Xnp1, Ynp1 = get_phy_char(npatchs, advec, which_advec, tstep, dt)
        # Computing errors:
        list_errs = []
        comp_err_time(Xnp1, Ynp1, listZ, list_errs, tval = tstep)
        list_err_inf.append(list_errs[0])
        list_err_l2.append(list_errs[1])
        list_minval.append(list_errs[2])
        list_maxval.append(list_errs[3])
        list_mass.append(list_errs[4])
        # Printing some results
        print '= =========== TSTEP = ', tstep, "/", nstep, ' =========== ='
        maxerr = np.max(list_errs[0])
        npat_maxerr = list_errs[0].index(np.max(list_errs[0]))
        print " --> For npat =", npat_maxerr, " maximum err l_inf =", maxerr


# -------------------------------------------
# Computing exact solution and plotting error :
Xnp1, Ynp1 = get_phy_char(npatchs, advec, which_advec, tstep, dt)
comp_err_time(Xnp1, Ynp1, \
              listZ, list_errs,
              show = True, plot = True, block=True,
              tval = tstep)
plot_nrb_dens(zn, show = True, save = False)

plot_errors([list_err_inf, list_err_l2, list_minval, list_maxval, list_mass])
