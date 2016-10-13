"""
Module to compute characteristics' origin both in physical
and logical domain.
"""

from globals_variables import *
from geometry import geo
from geometry import npatchs
from geometry import jac
from geometry import X_mat, Y_mat
import connectivity as conn


def get_pat_char(eta1, eta2, advec, dt):
    """
    Computes the characteristics in logical domain.

    Note:
        To compute characteristic's feet we have to
        solve the following ODE :
        d(Eta)/dt = 1/sqrtg rot(psi)
        with i.c. Eta(s) = eta(0)
        where
        sqrtg = det(J) = det (Jacobian)
        rot(Psi) = J^-1 * A * sqrtg
    Args:
        eta1: vector containing the coordinates of eta1
        eta2: vector containing the coordinates of eta2
        advec: advection coefficient for every direction, patch and point
        dt: time step

    Returns:
        char_eta1: 1st coordinate of origin of characteristc in logical domain
        char_eta2: 2nd coordinate of origin of characteristc in logical domain
    """

    char_eta1 = np.zeros((npatchs, NPTS1, NPTS2))
    char_eta2 = np.zeros((npatchs, NPTS1, NPTS2))
    where_char = np.zeros((npatchs, NPTS1, NPTS2), dtype=np.int)

    for npat, nrb in enumerate(geo):
        print ""
        print "For patch ", npat, " :"
        eta1_mat, eta2_mat = my_meshgrid(eta1, eta2)
        # Definition of right hand side of ODE :
        rhs = lambda xatn, yatn: derivs_eta(geo, npat,
                                            xatn, yatn,
                                            [advec[0][npat], advec[1][npat]])
        # We solve the ode using a Runge Kutta method of order 4
        char_eta1[npat], char_eta2[npat], where_char[npat] \
            = rk4(npat, eta1_mat, eta2_mat, dt, rhs)
        print "   * char_eta1 max et min =", \
                np.max(char_eta1[npat]), np.min(char_eta1[npat])
        print "   * char_eta2 max et min =", \
                np.max(char_eta2[npat]), np.min(char_eta2[npat])

    return char_eta1, char_eta2, where_char


def rk4(npat, xn, yn, dt, rhs_func):
    """
    Runge Kutta method of order 4 applied to ODE for
    characteristics (see get_pat_char).
    This function returns the results constrained to the patch.
    No out of domain particle.

    Args:
        npat: patch index
        xn: 1st coordinates at time tn
        yn: 2nd coordinates at time tn
        dt: time step
        rhs_func: function given the RHS of the ODE

    Returns:
        xnp1, ynp1: coordinates of characteristics origin
        advec_done_x, advec_done_y: advection done (=rhs unless
                                    particle out of domain)
    """

    # Runge-Kutta Method of order 4 :
    k1x, k1y = rhs_func(xn, yn)

    advec_x = 0.5 * dt * k1x
    advec_y = 0.5 * dt * k1y
    xn_05 = xn - advec_x
    yn_05 = yn - advec_y
    last_advec_percent = np.zeros_like(xn)
    where_char = np.zeros_like(xn, dtype=np.int) + npat
    contain_particles(xn, yn,
                      [advec_x, advec_y],
                      xn_05, yn_05,
                      where_char, last_advec_percent)

    k2x, k2y = rhs_func(xn_05, yn_05)

    advec_x = 0.5 * dt * k2x
    advec_y = 0.5 * dt * k2y
    xn_15 = xn - advec_x
    yn_15 = yn - advec_y
    last_advec_percent = np.zeros_like(xn)
    where_char = np.zeros_like(xn, dtype=np.int) + npat
    contain_particles(xn, yn,
                      [advec_x, advec_y],
                      xn_15, yn_15,
                      where_char, last_advec_percent)

    k3x, k3y = rhs_func(xn_15, yn_15)

    advec_x = dt * k3x
    advec_y = dt * k3y
    xn_25 = xn - advec_x
    yn_25 = yn - advec_y
    last_advec_percent = np.zeros_like(xn)
    where_char = np.zeros_like(xn, dtype=np.int) + npat
    contain_particles(xn, yn,
                      [advec_x, advec_y],
                      xn_25, yn_25,
                      where_char, last_advec_percent)

    k4x, k4y = rhs_func(xn_25, yn_25)

    advec_x = 0.5 * dt / 3. * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
    advec_y = 0.5 * dt / 3. * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)

    xnp1 = xn - advec_x
    ynp1 = yn - advec_y
    last_advec_percent = np.zeros_like(xn)
    where_char = np.zeros_like(xn, dtype=np.int) + npat
    if DEBUG_MODE:
        print "...................... last call ....................."
    contain_particles(xn, yn,
                      [advec_x, advec_y],
                      xnp1, ynp1,
                      where_char, last_advec_percent)

    print " last advec percent = ", last_advec_percent[np.where(last_advec_percent < 1)]
    print "where = ", np.where(last_advec_percent < 1)
    return xnp1, ynp1, where_char


def contain_particles(eta1_mat, eta2_mat,
                      advec_coeffs,
                      eta1_orig, eta2_orig,
                      where_char, last_advec_percent):
    """
    Resets values of characteristics origin such that their
    coordinates remain on the patch.
    Furthermore it saves the value of the advection done in eta1 and eta2.
    RECURSIVE FUNCTION.

    Notes:
        rhs_eta1, rhs_eta2, and dt, are only needed to compute
        the advection done (when not the characteristics origin
        leaves the patch). The value of the advection done in
        each direction is computed and stored in
        advec_done_1 and advec_done_2.
    Args:
        eta1_mat: 1st coordinates on patch (of all points)
        eta2_mat: 2nd coordinates on patch (of all points)
        rhs_eta1: 1st coordinate of RHS of the ODE to be solved
                  for characteristics origin
        rhs_eta2: 2nd coordinate of RHS of the ODE to be solved
                  for characteristics origin
        advec_done_1: [INOUT] value of the advection done in eta1
        advec_done_2: [INOUT] value of the advection done in eta2
        eta1_orig: [INOUT] 1st coordinate of characteristic's origin
                   truncated so that it stays in patch
        eta2_orig: [INOUT] 2nd coordinate of characteristic's origin
                   truncated so that it stays in patch
        where_char: [INOUT] index of patch where the origin of the
                    characteristic is situated.
    Returns:
    """

    # Separating advection:
    [advec_coef1, advec_coef2] = advec_coeffs

    # Rounding off small values to avoid useless loops :
    eta1_orig[np.where(abs(eta1_orig) < epsilon2)] = 0.
    eta2_orig[np.where(abs(eta2_orig) < epsilon2)] = 0.
    eta1_orig[np.where(abs(1. - eta1_orig) < epsilon2)]=1.
    eta2_orig[np.where(abs(1. - eta2_orig) < epsilon2)]=1.
    advec_coef1[np.where(abs(advec_coef1) < epsilon2)] = 0.
    advec_coef2[np.where(abs(advec_coef2) < epsilon2)] = 0.
    #...................................................

    if DEBUG_MODE :
        print "___________________________________________"
        print "----- eta mat"
        print eta1_mat[0, 48]
        print eta2_mat[0, 48]
        print "---------- advec"
        print advec_coef1[0, 48]
        print advec_coef2[0, 48]
        print "---------- origins"
        print eta1_orig[0, 48]
        print eta2_orig[0, 48]
        print "---------- where"
        print where_char[0, 48]
        print "---------- last percent"
        print last_advec_percent[0, 48]
        print "___________________________________________"


    # For particles that stay in the domain
    in_indices = np.where((eta1_orig >= 0.) & (eta1_orig <= 1.)
                          & (eta2_orig >= 0.) & (eta2_orig <= 1.))

    # Nothing to for particles that are in domain
    last_advec_percent[in_indices] = 1
    if (last_advec_percent == 1).all() and DEBUG_MODE :
        print "results ========="
        print "origins of chars: "
        print eta1_orig
        print eta2_orig
        print "in the pats:", where_char
        print "percent :", last_advec_percent
        print "===== leaving 1 ========"
        return

    # For particles that go out on the x axis, below 0 :
    where_out = np.where((eta1_orig < 0.) & (last_advec_percent < 1.))
    if np.size(where_out) > 0:
        contain_particles_1D(where_char, where_out, last_advec_percent,\
                             False, True, \
                             eta1_mat, eta2_mat, \
                             advec_coef1, advec_coef2, \
                             eta1_orig, eta2_orig, \
                             1)

    # For particles that go out on the x axis, above 1 :
    where_out = np.where((eta1_orig - 1. > 0.) & (last_advec_percent < 1.))
    if np.size(where_out) > 0:
        contain_particles_1D(where_char, where_out, last_advec_percent,\
                             True, True, \
                             eta1_mat, eta2_mat, \
                             advec_coef1, advec_coef2, \
                             eta1_orig, eta2_orig, \
                             3)

    # For particles that go out on the y axis, below 0 :
    where_out = np.where((eta2_orig < 0.) & (last_advec_percent < 1.))
    if np.size(where_out) > 0:
        contain_particles_1D(where_char, where_out, last_advec_percent,\
                             False, False, \
                             eta2_mat, eta1_mat, \
                             advec_coef2, advec_coef1, \
                             eta2_orig, eta1_orig, \
                             0)

    # For particles that go out on the y axis, above 1 :
    where_out = np.where((eta2_orig - 1.0 > 0.) & (last_advec_percent < 1.))
    if np.size(where_out) > 0:
        contain_particles_1D(where_char, where_out, last_advec_percent,\
                             True, False, \
                             eta2_mat, eta1_mat, \
                             advec_coef2, advec_coef1, \
                             eta2_orig, eta1_orig, \
                             2)

    eta1_orig[np.where(abs(eta1_orig) < epsilon2)] = 0.
    eta2_orig[np.where(abs(eta2_orig) < epsilon2)] = 0.

    if DEBUG_MODE:
        print "results ========="
        print "origins of chars: "
        print eta1_orig
        print eta2_orig
        print "in the pats:", where_char
        print "percent :", last_advec_percent
        print "________leaving 2__________"

    out_indices = np.where((eta1_orig < 0.) | (eta1_orig > 1.)
                        | (eta2_orig < 0.) | (eta2_orig > 1.))

    if (np.size(out_indices) > 0):
        print "ERROR in contain_particles(): "\
            +"some particles are still outside domain!"
        print "indices = ", out_indices
        print eta1_orig[out_indices]
        print eta2_orig[out_indices]
        print ""
        import os, sys
        sys.exit("STOP")

def contain_particles_1D(where_char, where_out, last_advec_percent,\
                         is_above1, is_eta1,\
                         eta_out, eta_in, \
                         advec_out, advec_in,\
                         eta_out_orig, eta_in_orig, \
                         face):

    # We compute the actual percent we could advect from
    current_percent = (eta_out[where_out] - 1.*is_above1)/ advec_out[where_out]
    # Checking if the particle is not still out of the patch
    temp = eta_in[where_out] - current_percent * advec_in[where_out]
    where_out2 = where_out[0][np.where((temp >= 0.) & (temp <=1.))[0]]
    if np.shape(where_out)[0]!=1:
        where_out3 = where_out[1][np.where((temp >= 0.) & (temp <=1.))[0]]
        where_out = (where_out2, where_out3)
    else:
        where_out = (where_out2,)
    # updating the percent:
    current_percent = (eta_out[where_out] -1.*is_above1) / advec_out[where_out]
    # ....
    last_advec_percent[where_out] += current_percent

    # We compute where the particle stayed
    eta_out_orig[where_out] = 1.*is_above1
    eta_in_orig[where_out] = eta_in[where_out] \
                           - current_percent * advec_in[where_out]

    # Looking for the neigh. patch and transforming the coord to that patch
    [list_pats, list_faces] = conn.connectivity(where_char[where_out], face)

    # We treat external external boundary conditions here,
    # ie when [npat, face] = [npat', face']
    where_not_dirichlet = np.where((where_char[where_out] != list_pats)
                                   & (np.asarray(list_faces) != face))[0]
    where_is_dirichlet = np.where((where_char[where_out] == list_pats)
                                  & (np.asarray(list_faces) == face))[0]
    if np.shape(where_out)[0]!=1:
        where_dirichlet = [where_out[0][where_is_dirichlet],
                           where_out[1][where_is_dirichlet]]
        where_out = [where_out[0][where_not_dirichlet],
                     where_out[1][where_not_dirichlet]]
    else:
        where_dirichlet = where_out[0][where_is_dirichlet]
        np.reshape(where_dirichlet,(1, np.size(where_dirichlet)))
        where_out = where_out[0][where_not_dirichlet]
        np.reshape(where_out,(1, np.size(where_out)))

    # Treating DIRICHLET BC .........................
    # We can not advec any further so we impose the percent to complete state
    last_advec_percent[where_dirichlet] = 1.
    #................................................

    if np.size(where_out) > 0:
        where_orig = where_char[where_out]
        list_pats2 = [list_pats[val0] for index0, val0 \
                      in enumerate(where_not_dirichlet)]
        where_char[where_out] = list_pats2

        [eta1_out_xmin, eta2_out_xmin] = \
            conn.transform_patch_coordinates(eta_in_orig[where_out], list_faces)

        if (is_eta1) :
            [advec1, advec2] = \
                    conn.transform_advection(\
                    advec_out[where_out], advec_in[where_out],
                    where_orig, where_char[where_out],
                    eta_out_orig[where_out], eta_in_orig[where_out], \
                    eta1_out_xmin, eta2_out_xmin)
            # We get the origin point
            eta_out_orig[where_out] = eta1_out_xmin
            eta_in_orig[where_out] = eta2_out_xmin
        else :
            [advec1, advec2] = \
                    conn.transform_advection(\
                    advec_in[where_out], advec_out[where_out],
                    where_orig, where_char[where_out],
                    eta_in_orig[where_out], eta_out_orig[where_out], \
                    eta1_out_xmin, eta2_out_xmin)
            # We get the origin point
            eta_in_orig[where_out] = eta1_out_xmin
            eta_out_orig[where_out] = eta2_out_xmin

        # Now we advect with the remaining percentage left
        eta1_out_xmin += -advec1 * (1. - last_advec_percent[where_out])
        eta2_out_xmin += -advec2 * (1. - last_advec_percent[where_out])
        # and we contain the particles again
        this_percent = last_advec_percent[where_out]
        if (this_percent < 0.).any():
            import sys
            sys.exit("In contain_particles_1D: negative percentage found")

        if (is_eta1):
            e1o = eta_out_orig[where_out]
            e2o = eta_in_orig[where_out]
            wc = where_char[where_out]

            contain_particles(e1o, e2o, \
                          [advec1, advec2], \
                          eta1_out_xmin, eta2_out_xmin, \
                          wc, this_percent)

            # We can now replace with new values:
            eta_out_orig[where_out] = eta1_out_xmin
            eta_in_orig[where_out] = eta2_out_xmin
            where_char[where_out] = wc
            last_advec_percent[where_out] = this_percent

        else:
            e1o = eta_in_orig[where_out]
            e2o = eta_out_orig[where_out]
            wc = where_char[where_out]

            contain_particles(e1o, e2o, \
                              [advec1, advec2], \
                              eta1_out_xmin, eta2_out_xmin, \
                              wc, this_percent)
            # We can now replace with new values:
            eta_in_orig[where_out] = eta1_out_xmin
            eta_out_orig[where_out] = eta2_out_xmin
            where_char[where_out] = wc
            last_advec_percent[where_out] = this_percent
    # else:
    #     print ""
    #     print " where out is empty"
    #     print ""



def derivs_eta(geo, npat, eta1_mat, eta2_mat, advec):

    # Getting the jacobian
    [d1F1, d2F1, d1F2, d2F2] = jacobian(geo, npat, eta1_mat, eta2_mat)

    # Computing the jacobian determinant
    sqrt_g = d1F1 * d2F2 - d1F2 * d2F1
    # print "++++++ jacobien = ", sqrt_g

    # Approximating the 0 values of sqrtg to avoid Nans
    # if (np.size(np.where(sqrt_g == 0.)) > 0) :
    #     print "heeeere", np.where(sqrt_g == 0.)
    #     sqrt_g[0,0] = 0.5*(sqrt_g[1,0] + sqrt_g[0,1])
    #     sqrt_g[NPTS1-1,0] = 0.5*(sqrt_g[NPTS1-2,0] + sqrt_g[NPTS1-1,1])
    #     sqrt_g[0,NPTS2-1] = 0.5*(sqrt_g[1,NPTS2-1] + sqrt_g[0,NPTS2-2])
    #     sqrt_g[0,NPTS1-1] = 0.5*(sqrt_g[NPTS1-2,NPTS2-1] + sqrt_g[NPTS1-1,NPTS2-2])
    # print np.where(sqrt_g == 0.)

    wjm = np.where(abs(sqrt_g) < epsilon2)
    if (np.size(wjm) > 0):
        print "Warning : singular point"
        print np.size(wjm)
        # import os, sys
        # sys.exit("STOP")
        sqrt_g[wjm] = epsilon2

    # Calculating the value of the second part of the MoC
    rhs1 = (advec[0] * d2F2 - advec[1] * d2F1) / sqrt_g
    rhs2 = (advec[1] * d1F1 - advec[0] * d1F2) / sqrt_g

    # print "a 0 =", advec[0]
    # print "d2f2 =", d2F2
    # print "a 1 =", advec[1]
    # print "d2f1 =", d2F1
    # print "sqrtg = ", sqrt_g
    # print "advec1 in pat ",rhs1
    # print "advec2 in pat ",rhs2
    return rhs1, rhs2


def jacobian(geo, npat, eta1_mat, eta2_mat):
    """
    Computes jacobian in points given.

    Args:
        nrb: Contains info about transformations and patches, given by igakit
        eta1_mat: matrix containing eta1 coordinates
        eta2_mat: matrix containing eta2 coordinates

    Returns:
        [d1F1, d2F1, d1F2, d2F2]: list containing the values
        of the jacobian matrix at given points.
    """
    u = eta1_mat[:, 0]
    v = eta2_mat[0, :]

    wum = np.where(u < 0.)
    wup = np.where(u > 1.)
    if ((np.size(wum) > 0) or (np.size(wup) > 0)):
        print "Warning: in jacobian() from charac_feet.py:"
        print "         Found a value outside [0,1] in vector u"
        u[wum] = 0.
        u[wup] = 1.
    wvm = np.where(v < 0.)
    wvp = np.where(v > 1.)
    if ((np.size(wvm) > 0) or (np.size(wvp) > 0)):
        print "Warning: in jacobian() from charac_feet.py:"
        print "         Found a value outside [0,1] in vector v"
        v[wvm] = 0.
        v[wvp] = 1.

    d1F1 = np.zeros((NPTS1, NPTS2))
    d2F1 = np.zeros((NPTS1, NPTS2))
    d1F2 = np.zeros((NPTS1, NPTS2))
    d2F2 = np.zeros((NPTS1, NPTS2))

    # Getting the derivatives
    d1F1 = jac[npat,0,:,:]
    d2F1 = jac[npat,1,:,:]
    d1F2 = jac[npat,2,:,:]
    d2F2 = jac[npat,3,:,:]

    # Getting rid of close to 0 values
    d1F1[np.where(abs(d1F1) <= small_epsilon)] = 0.0
    d2F1[np.where(abs(d2F1) <= small_epsilon)] = 0.0
    d1F2[np.where(abs(d1F2) <= small_epsilon)] = 0.0
    d2F2[np.where(abs(d2F2) <= small_epsilon)] = 0.0

    return [d1F1, d2F1, d1F2, d2F2]


def get_phy_char(npatchs, advec, which_advec, tstep, dt):
    """
    Computes characteristics' origin in the physical domain.

    Args:
        npatchs: number of patches
        advec: advection coefficients on each coordinate of each patch
        which_advec: type of advection (int)
        tstep: number of iteration being done in time
        dt: step in time

    Returns:
        xnp1: 1st coordinates of characteristics origin at time tstep
        ynp1: 2nd coordinates of characteristics origin at time tstep
    """
    xnp1 = np.zeros((npatchs, NPTS1, NPTS2))
    ynp1 = np.zeros((npatchs, NPTS1, NPTS2))

    for npat in range(npatchs):
        if which_advec == 0:
            xnp1[npat] = X_mat[npat] - advec[0][npat] * tstep * dt
            ynp1[npat] = Y_mat[npat] - advec[1][npat] * tstep * dt

        elif which_advec == 1:
            r = np.sqrt(X_mat[npat] ** 2 + Y_mat[npat] ** 2)
            th = np.arctan2(Y_mat[npat], X_mat[npat])
            # TODO : maybe it's better to use directly X_mat and Y_mat
            xnp1[npat] = r * np.cos(-2. * np.pi * tstep * dt + th)
            ynp1[npat] = r * np.sin(-2. * np.pi * tstep * dt + th)

        elif which_advec == 2:
            r = np.sqrt(advec[0][npat] ** 2 + advec[1][npat] ** 2)
            th = np.arctan2(advec[0][npat], -advec[1][npat])
            # TODO : maybe it's better to use directly X_mat and Y_mat
            xnp1[npat] = r * np.sin((tstep) * dt + np.pi / 2. - th) + centera
            ynp1[npat] = r * np.cos((tstep) * dt + np.pi / 2. - th) + centerb

        else:
            print "ERROR: in get_phy_char() not known advection."
            import os, sys
            sys.exit("STOP")

    return xnp1, ynp1
