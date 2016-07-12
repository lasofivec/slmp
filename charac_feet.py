from globals_variables import *
import interpol as inter


# ============================================================
#
# Module to compute characteristics' origin both in physical
# and logical domain.
#
# ============================================================



def get_phy_char(X_mat, Y_mat, npatchs, advec, which_advec, tstep, dt):
    """
    Computes characteristics' origin in the physical domain.

    Args:
        X_mat: 1st coordinates of mesh points in physical domain
        Y_mat: 2nd coordinates of mesh points in physical domain
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


def get_pat_char(geo, eta1, eta2, advec, dt):
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
        geo: Set of all the patches, given by igakit
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
        rhs = lambda xatn, yatn: derivs_eta(geo, npat, xatn, yatn, [advec[0][npat], advec[1][npat]])
        # we solve the ode using a Runge Kutta method of order 4
        char_eta1[npat], char_eta2[npat], where_char[npat] = rk4(npat, eta1_mat, eta2_mat, dt, rhs)
        print "   * char_eta1 max et min =", np.max(char_eta1[npat]), np.min(char_eta1[npat])
        print "   * char_eta2 max et min =", np.max(char_eta2[npat]), np.min(char_eta2[npat])


    return char_eta1, char_eta2, where_char


def rk4(npat, xn, yn, dt, rhs_func):
    """
    Runge Kutta method of order 4 applied to ODE for characteristics (see get_pat_char)
    This function returns the results constrained to the patch. No out of domain particle.

    Args:
        npat: patch index
        xn: 1st coordinates at time tn
        yn: 2nd coordinates at time tn
        dt: time step
        rhs_func: function given the RHS of the ODE

    Returns:

        xnp1, ynp1: coordinates of characteristics origin
        advec_done_x, advec_done_y: advection done (=rhs unless particle out of domain)
    """

    # Runge-Kutta Method of order 4 :
    k1x, k1y = rhs_func(xn, yn)

    advec_x = 0.5 * dt * k1x
    advec_y = 0.5 * dt * k1y
    xn_05 = xn - advec_x
    yn_05 = yn - advec_y
    last_advec_percent = np.zeros_like(xn)
    where_char = np.zeros_like(xn, dtype=np.int) + npat
    contain_particles(xn, yn, [advec_x, advec_y], xn_05, yn_05, where_char, last_advec_percent)
    k2x, k2y = rhs_func(xn_05, yn_05)

    advec_x = 0.5 * dt * k2x
    advec_y = 0.5 * dt * k2y
    xn_15 = xn - advec_x
    yn_15 = yn - advec_y
    last_advec_percent = np.zeros_like(xn)
    where_char = np.zeros_like(xn, dtype=np.int) + npat
    contain_particles(xn, yn, [advec_x, advec_y], xn_15, yn_15, where_char, last_advec_percent)
    k3x, k3y = rhs_func(xn_15, yn_15)

    advec_x = dt * k3x
    advec_y = dt * k3y
    xn_25 = xn - advec_x
    yn_25 = yn - advec_y
    last_advec_percent = np.zeros_like(xn)
    where_char = np.zeros_like(xn, dtype=np.int) + npat
    contain_particles(xn, yn, [advec_x, advec_y], xn_25, yn_25, where_char, last_advec_percent)
    k4x, k4y = rhs_func(xn_25, yn_25)

    advec_x = 0.5 * dt / 3. * (k1x + 2.0 * k2x + 2.0 * k3x + k4x)
    advec_y = 0.5 * dt / 3. * (k1y + 2.0 * k2y + 2.0 * k3y + k4y)
    print ",,,,,,,,,,,,,,, last call ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"
    xnp1 = xn - advec_x
    ynp1 = yn - advec_y
    last_advec_percent = np.zeros_like(xn)
    where_char = np.zeros_like(xn, dtype=np.int) + npat
    contain_particles(xn, yn, [advec_x, advec_y], xnp1, ynp1, where_char, last_advec_percent)
    print ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"
    return xnp1, ynp1, where_char


def contain_particles(eta1_mat, eta2_mat, advec_coeffs, eta1_orig, eta2_orig, where_char, last_advec_percent):
    """
    Resets values of characteristics origin such that their coordinates remain on the patch.
    Furthermore it saves the value of the advection done in eta1 and eta2.
    RECURSIVE FUNCTION.

    Notes:
        rhs_eta1, rhs_eta2, and dt, are only need to compute the advection done (when not the characteristics origin
        leaves the patch). The value of the advection done in each direction is computed and stored
        in advec_done_1 and advec_done_2.
    Args:
        eta1_mat: 1st coordinates on patch (of all points)
        eta2_mat: 2nd coordinates on patch (of all points)
        rhs_eta1: 1st coordinate of RHS of the ODE to be solved for characteristics origin
        rhs_eta2: 2nd coordinate of RHS of the ODE to be solved for characteristics origin
        advec_done_1: [INOUT] value of the advection done in eta1
        advec_done_2: [INOUT] value of the advection done in eta2
        eta1_orig: [INOUT] 1st coordinate of characteristic's origin truncated so that it stays in patch
        eta2_orig: [INOUT] 2nd coordinate of characteristic's origin truncated so that it stays in patch
        where_char: [INOUT] index of patch where the origin of the characteristic is situated.

    Returns:
    """

    # Separating advection:
    [advec_coef1, advec_coef2] = advec_coeffs

    # # Rounding results:
    eta1_orig[np.where(abs(eta1_orig) < epsilon2)]=0.
    eta2_orig[np.where(abs(eta2_orig) < epsilon2)]=0.
    eta1_orig[np.where(abs(1. - eta1_orig) < epsilon2)]=1.
    eta2_orig[np.where(abs(1. - eta2_orig) < epsilon2)]=1.

    if DEBUG_MODE:
        print "---------- eta mat"
        print eta1_mat
        print eta2_mat
        print "---------- advec"
        print advec_coef1
        print advec_coef2
        print "---------- origins"
        print eta1_orig
        print eta2_orig
        print "---------- where"
        print where_char
        print "---------- last percent"
        print last_advec_percent
        
    # For particles that stay in the domain
    in_indices = np.where((eta1_orig >= 0.) & (eta1_orig <= 1.) & (eta2_orig >= 0.) & (eta2_orig <= 1.))
    # Nothing to for particles that are in domain
    last_advec_percent[in_indices] = 100
    if (last_advec_percent == 100).all() :
        return

    # For particles that go out on the x axis, below 0 :
    where_out = np.where(eta1_orig < 0.)
    if np.size(where_out) > 0:
        # We compute the actual percent we could advect from
        current_percent = eta1_mat[where_out] / advec_coef1[where_out]
        # Checking if the particle is not still out of the patch
        temp = eta2_mat[where_out] - current_percent * advec_coef2[where_out]
        where_out2 = where_out[0][np.where((temp >= 0.) & (temp <=1.))[0]]
        if np.shape(where_out)[0]!=1:
            where_out3 = where_out[1][np.where((temp >= 0.) & (temp <=1.))[0]]
            where_out = (where_out2, where_out3)
        else:
            where_out = (where_out2,)
        current_percent = eta1_mat[where_out] / advec_coef1[where_out]
        # ....
        last_advec_percent[where_out] += current_percent

        # We compute where the particle stayed
        eta1_orig[where_out] = 0.
        eta2_orig[where_out] = eta2_mat[where_out] - current_percent * advec_coef2[where_out]

        # Looking for the neigh. patch and transforming the coord to that patch
        [list_pats, list_faces]= inter.connectivity(where_char[where_out], 1)
        
        # We treat external external boundary conditions here, ie when [npat, face] = [npat', face']
        where_not_dirichlet = np.where((where_char[where_out] != list_pats) & (np.asarray(list_faces) != 1))[0]
        if np.shape(where_out)[0]!=1:
            where_out = [where_out[0][where_not_dirichlet], where_out[1][where_not_dirichlet]]
        else:
            where_out = where_out[0][where_not_dirichlet]
            np.reshape(where_out,(1, np.size(where_out)))

        if np.size(where_out) > 0:
            where_orig = where_char[where_out]
            list_pats2 = [list_pats[val0] for index0, val0 in enumerate(where_not_dirichlet)]
            where_char[where_out] = list_pats2

            [eta1_out_xmin, eta2_out_xmin] = \
                inter.transform_patch_coordinates(eta2_orig[where_out], list_faces)
            
            [advec1, advec2] = \
                inter.transform_advection(advec_coef1[where_out], advec_coef2[where_out], where_orig, where_char[where_out], eta1_orig[where_out], eta2_orig[where_out], eta1_out_xmin, eta2_out_xmin)
            
            # We get the origin point
            eta1_orig[where_out] = eta1_out_xmin
            eta2_orig[where_out] = eta2_out_xmin
            # Now we advect with the remaining percentage left
            eta1_out_xmin += -advec1 * (1. - last_advec_percent[where_out])
            eta2_out_xmin += -advec2 * (1. - last_advec_percent[where_out])
            # and we contain the particles again
            this_percent = last_advec_percent[where_out]
            if (this_percent < 0.).any():
                print "mats  = ", eta1_mat[where_out], eta2_mat[where_out]
                print "orig  = ", eta1_orig[where_out], eta2_orig[where_out]
                print "advec = ",advec1, advec2
                print "chara = ", eta1_out_xmin, eta2_out_xmin
                print "this percent", this_percent
                import sys
                sys.exit("In eta1 < 0: negative percentage found")
            contain_particles(eta1_orig[where_out], eta2_orig[where_out], \
                              [advec1, advec2], \
                              eta1_out_xmin, eta2_out_xmin, \
                              where_char[where_out], this_percent)
            # We can now replace with new values:
            eta1_orig[where_out] = eta1_out_xmin
            eta2_orig[where_out] = eta2_out_xmin
            last_advec_percent[where_out] = this_percent
            if DEBUG_MODE:
                print "results ========="
                print "where out =", where_out
                print "origins of chars: ", eta1_orig, eta2_orig
                print "in the pats:", where_char
                print "percent :", this_percent
                print "advec =", advec1

    # For particles that go out on the x axis, above 1 :
    where_out = np.where(eta1_orig - 1. > 0.)
    if np.size(where_out) > 0:
        # We compute the actual percent we could advect from
        current_percent = (eta1_mat[where_out] - 1.) / advec_coef1[where_out]
        # Checking if the particle is not still out of the patch
        temp = eta2_mat[where_out] - current_percent * advec_coef2[where_out]
        where_out2 = where_out[0][np.where((temp >= 0.) & (temp <=1.))[0]]
        if np.shape(where_out)[0]!=1:
            where_out3 = where_out[1][np.where((temp >= 0.) & (temp <=1.))[0]]
            where_out = (where_out2, where_out3)
        else:
            where_out = (where_out2,)
        current_percent = (eta1_mat[where_out] - 1.) / advec_coef1[where_out]
        # ....
        last_advec_percent[where_out] += current_percent

        # We compute where the particle stayed
        eta1_orig[where_out] = 1.
        eta2_orig[where_out] = eta2_mat[where_out] - current_percent * advec_coef2[where_out]

        # Looking for the neighbouring patch and transforming the coordinates to that patch
        [list_pats, list_faces]= inter.connectivity(where_char[where_out], 3)

        # We treat external external boundary conditions here, ie when [npat, face] = [npat', face']
        where_not_dirichlet = np.where((where_char[where_out] != list_pats) & (np.asarray(list_faces) != 3))[0]
        if np.shape(where_out)[0]!=1:
            where_out = [where_out[0][where_not_dirichlet], where_out[1][where_not_dirichlet]]
        else:
            where_out = where_out[0][where_not_dirichlet]
            np.reshape(where_out,(1, np.size(where_out)))

        if np.size(where_out) > 0:
            where_orig = where_char[where_out]
            list_pats2 = [list_pats[val0] for index0, val0 in enumerate(where_not_dirichlet)]
            where_char[where_out] = list_pats
            [eta1_out_xmin, eta2_out_xmin] = \
                inter.transform_patch_coordinates(eta2_orig[where_out], list_faces)

            [advec1, advec2] = \
                inter.transform_advection(advec_coef1[where_out], advec_coef2[where_out], where_orig, where_char[where_out], eta1_orig[where_out], eta2_orig[where_out], eta1_out_xmin, eta2_out_xmin)

            # We get the origin point
            eta1_orig[where_out] = eta1_out_xmin
            eta2_orig[where_out] = eta2_out_xmin
            # Now we advect with the remaining percentage left
            eta1_out_xmin += - advec1 * (1. - last_advec_percent[where_out])
            eta2_out_xmin += - advec2 * (1. - last_advec_percent[where_out])
            # and we contain the particles again
            this_percent = last_advec_percent[where_out]
            if (this_percent < 0.).any():
                print ".................. "
                print "mat  = ", eta1_mat[where_out], eta2_mat[where_out]
                print "orig  = ", eta1_orig[where_out], eta2_orig[where_out]
                print "advec = ",advec1, advec2
                print "chara = ", eta1_out_xmin, eta2_out_xmin
                print "this percent", this_percent
                print "current percent", current_percent
                import sys
                sys.exit("In eta1 - 1 > 0: negative percentage found")
            contain_particles(eta1_orig[where_out], eta2_orig[where_out], \
                              [advec1, advec2], \
                              eta1_out_xmin, eta2_out_xmin, \
                              where_char[where_out], this_percent)
            # We can now replace with new values:
            eta1_orig[where_out] = eta1_out_xmin
            eta2_orig[where_out] = eta2_out_xmin
            last_advec_percent[where_out] = this_percent

    # For particles that go out on the y axis, below 0 :
    where_out = np.where(eta2_orig < 0.)
    if np.size(where_out) > 0:
        # We compute the actual percent we could advect from
        current_percent = eta2_mat[where_out] / advec_coef2[where_out]
        last_advec_percent[where_out] += current_percent
        # We compute where the particle stayed
        eta1_orig[where_out] = eta1_mat[where_out] - current_percent * advec_coef1[where_out]
        eta2_orig[where_out] = 0.
        
        # Looking for the neighbouring patch and transforming the coordinates to that patch
        [list_pats, list_faces] = inter.connectivity(where_char[where_out], 0)

        # We treat external external boundary conditions here, ie when [npat, face] = [npat', face']
        where_not_dirichlet = np.where((where_char[where_out] != list_pats) & (np.asarray(list_faces) != 0))[0]
        if np.shape(where_out)[0]!=1:
            where_out = [where_out[0][where_not_dirichlet], where_out[1][where_not_dirichlet]]
        else:
            where_out = where_out[0][where_not_dirichlet]
            np.reshape(where_out,(1, np.size(where_out)))
                
        if np.size(where_out) > 0:
            where_orig = where_char[where_out]
            list_pats2 = [list_pats[val0] for index0, val0 in enumerate(where_not_dirichlet)]
            where_char[where_out] = list_pats2

            [eta1_out_xmin, eta2_out_xmin] = \
                inter.transform_patch_coordinates(eta1_orig[where_out], list_faces)

            [advec1, advec2] = \
                inter.transform_advection(advec_coef1[where_out], advec_coef2[where_out], where_orig, where_char[where_out], eta1_orig[where_out], eta2_orig[where_out], eta1_out_xmin, eta2_out_xmin)

            # We get the origin point
            eta1_orig[where_out] = eta1_out_xmin
            eta2_orig[where_out] = eta2_out_xmin
            # Now we advect with the remaining percentage left
            eta1_out_xmin += - advec1 * (1. - last_advec_percent[where_out])
            eta2_out_xmin += - advec2 * (1. - last_advec_percent[where_out])
            # and we contain the particles again
            this_percent = last_advec_percent[where_out]
            if (this_percent < 0.).any():
                print "orig  = ", eta1_orig[where_out], eta2_orig[where_out]
                print "advec = ",advec1, advec2
                print "chara = ", eta1_out_xmin, eta2_out_xmin
                print "this percent", this_percent
                import sys
                sys.exit("In eta2 < 0: negative percentage found")
            contain_particles(eta1_orig[where_out], eta2_orig[where_out], \
                              [advec1, advec2], \
                              eta1_out_xmin, eta2_out_xmin, \
                              where_char[where_out], this_percent)
            # We can now replace with new values:
            eta1_orig[where_out] = eta1_out_xmin
            eta2_orig[where_out] = eta2_out_xmin
            last_advec_percent[where_out] = this_percent

    # For particles that go out on the y axis, above 1 :
    where_out = np.where(eta2_orig - 1.0 > 0.)
    if np.size(where_out) > 0:
        # We compute the actual percent we could advect from
        current_percent = (eta2_mat[where_out] - 1.) / advec_coef2[where_out]
        last_advec_percent[where_out] += current_percent
        # We compute where the particle stayed
        eta1_orig[where_out] = eta1_mat[where_out] - current_percent * advec_coef1[where_out]
        eta2_orig[where_out] = 1.

        # Looking for the neighbouring patch and transforming the coordinates to that patch
        [list_pats, list_faces] = inter.connectivity(where_char[where_out], 2)

        # We treat external external boundary conditions here, ie when [npat, face] = [npat', face']
        where_not_dirichlet = np.where((where_char[where_out] != list_pats) & (np.asarray(list_faces) != 2))[0]
        if np.shape(where_out)[0]!=1:
            where_out = [where_out[0][where_not_dirichlet], where_out[1][where_not_dirichlet]]
        else:
            where_out = where_out[0][where_not_dirichlet]
            np.reshape(where_out,(1, np.size(where_out)))

        if np.size(where_out) > 0:
            where_orig = where_char[where_out]
            list_pats2 = [list_pats[val0] for index0, val0 in enumerate(where_not_dirichlet)]
            where_char[where_out] = list_pats2

            [eta1_out_xmin, eta2_out_xmin] = \
                inter.transform_patch_coordinates(eta1_orig[where_out], list_faces)
            
            [advec1, advec2] = \
                inter.transform_advection(advec_coef1[where_out], advec_coef2[where_out], where_orig, where_char[where_out], eta1_orig[where_out], eta2_orig[where_out], eta1_out_xmin, eta2_out_xmin)

            # We get the origin point
            eta1_orig[where_out] = eta1_out_xmin
            eta2_orig[where_out] = eta2_out_xmin
            # Now we advect with the remaining percentage left
            eta1_out_xmin += - advec1 * (1. - last_advec_percent[where_out])
            eta2_out_xmin += - advec2 * (1. - last_advec_percent[where_out])
            # and we contain the particles again
            this_percent = last_advec_percent[where_out]
            if (this_percent < 0.).any():
                print "orig  = ", eta1_orig[where_out], eta2_orig[where_out]
                print "advec = ",advec1, advec2
                print "chara = ", eta1_out_xmin, eta2_out_xmin
                print "this percent", this_percent
                import sys
                sys.exit("In eta2 - 1 > 0: negative percentage found")
            contain_particles(eta1_orig[where_out], eta2_orig[where_out], \
                              [advec1, advec2], \
                              eta1_out_xmin, eta2_out_xmin, \
                              where_char[where_out], this_percent)
            # We can now replace with new values:
            eta1_orig[where_out] = eta1_out_xmin
            eta2_orig[where_out] = eta2_out_xmin
            last_advec_percent[where_out] = this_percent

    out_indices = np.where((eta1_orig < 0.) & (eta1_orig > 1.) & (eta2_orig < 0.) & (eta2_orig > 1.))
    if (np.size(out_indices) > 0):
        print "ERROR in contain_particles(): some particles are still outside domain!"
        print "indices = ", out_indices
        print eta1_orig[out_indices]
        print eta2_orig[out_indices]
        print ""
        import os, sys
        sys.exit("STOP")


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

    if (np.size(np.where(sqrt_g == 0.)) > 0):
        print "Warning : singular point"
        print np.size(np.where(sqrt_g == 0. ))
        # import os, sys
        # sys.exit("STOP")
        sqrt_g[np.where(sqrt_g == 0.)] = epsilon

    # Calculating the value of the second part of the MoC
    rhs1 = (advec[0] * d2F2 - advec[1] * d2F1) / sqrt_g
    rhs2 = (advec[1] * d1F1 - advec[0] * d1F2) / sqrt_g

    return rhs1, rhs2


def jacobian(geo, npat, eta1_mat, eta2_mat):
    """
    Computes jacobian in points given.

    Args:
        nrb: Contains info about transformations and patches, given by igakit
        eta1_mat: matrix containing eta1 coordinates
        eta2_mat: matrix containing eta2 coordinates

    Returns:
        [d1F1, d2F1, d1F2, d2F2]: list containing the values of the jacobian matrix at given points.

    """
    u = eta1_mat[:, 0]
    v = eta2_mat[0, :]

    if ((np.size(np.where(u < 0.)) > 0) or (np.size(np.where(u > 1.)) > 0)):
        print "Warning: in jacobian() from charac_feet.py:"
        print "         Found a value outside [0,1] in vector u"
        u[np.where(u < 0.)] = 0.
        u[np.where(u > 1.)] = 1.
    if ((np.size(np.where(v < 0.)) > 0) or (np.size(np.where(v > 1.)) > 0)):
        print "Warning: in jacobian() from charac_feet.py:"
        print "         Found a value outside [0,1] in vector u"
        v[np.where(v < 0.)] = 0.
        v[np.where(v > 1.)] = 1.

    d1F1 = np.zeros((NPTS1, NPTS2))
    d2F1 = np.zeros((NPTS1, NPTS2))
    d1F2 = np.zeros((NPTS1, NPTS2))
    d2F2 = np.zeros((NPTS1, NPTS2))

    # Getting the derivatives
    if domain == DISK_5p_domain :
        jacobians = jacobian_function(npat, eta1_mat, eta2_mat)
        d1F1 = jacobians[0]
        d2F1 = jacobians[1]
        d1F2 = jacobians[2]
        d2F2 = jacobians[3]
    else :
        D = geo[npat].evaluate_deriv(u, v, nderiv=1)
        d1F1 = D[1, :, :, 0]
        d2F1 = D[2, :, :, 0]
        d1F2 = D[1, :, :, 1]
        d2F2 = D[2, :, :, 1]

    # Getting rid of close to 0 values
    d1F1[np.where(abs(d1F1) <= 10 ** -14)] = 0.0
    d2F1[np.where(abs(d2F1) <= 10 ** -14)] = 0.0
    d1F2[np.where(abs(d1F2) <= 10 ** -14)] = 0.0
    d2F2[np.where(abs(d2F2) <= 10 ** -14)] = 0.0

    return [d1F1, d2F1, d1F2, d2F2]


