def contain_particles_1D(where_out, \
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
    if np.shape(where_out)[0]!=1:
        where_out = [where_out[0][where_not_dirichlet],
                     where_out[1][where_not_dirichlet]]
    else:
        where_out = where_out[0][where_not_dirichlet]
        np.reshape(where_out,(1, np.size(where_out)))

    if np.size(where_out) > 0:
        where_orig = where_char[where_out]
        list_pats2 = [list_pats[val0] for index0, val0 \
                      in enumerate(where_not_dirichlet)]
        where_char[where_out] = list_pats2

        [eta1_out, eta2_out] = \
            conn.transform_patch_coordinates(eta_in_orig[where_out], list_faces)

        if (is_eta1) :
            [advec1, advec2] = \
                    conn.transform_advection(\
                    advec_out[where_out], advec_in[where_out],
                    where_orig, where_char[where_out],
                    eta_out_orig[where_out], eta_in_orig[where_out], \
                    eta1_out, eta2_out)
            # We get the origin point
            eta_out_orig[where_out] = eta1_out
            eta_in_orig[where_out] = eta2_out
        else :
            [advec1, advec2] = \
                    conn.transform_advection(\
                    advec_in[where_out], advec_out[where_out],
                    where_orig, where_char[where_out],
                    eta_in_orig[where_out], eta_out_orig[where_out], \
                    eta1_out, eta2_out)
            # We get the origin point
            eta_in_orig[where_out] = eta1_out
            eta_out_orig[where_out] = eta2_out

        # Now we advect with the remaining percentage left
        eta1_out += -advec1 * (1. - last_advec_percent[where_out])
        eta2_out += -advec2 * (1. - last_advec_percent[where_out])
        # and we contain the particles again
        this_percent = last_advec_percent[where_out]
        if (this_percent < 0.).any():
            import sys
            sys.exit("In contain_particles_1D: negative percentage found")

        if (is_eta1):
            contain_particles(eta_out_orig[where_out], eta_in_orig[where_out], \
                          [advec1, advec2], \
                          eta1_out, eta2_out, \
                          where_char[where_out], this_percent)
        else:
            contain_particles(eta_in_orig[where_out], eta_out_orig[where_out], \
                          [advec1, advec2], \
                          eta1_out, eta2_out, \
                          where_char[where_out], this_percent)

        # We can now replace with new values:
        eta1_orig[where_out] = eta1_out
        eta2_orig[where_out] = eta2_out
        last_advec_percent[where_out] = this_percent
