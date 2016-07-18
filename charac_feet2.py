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

    # For particles that stay in the domain
    in_indices = np.where((eta1_orig >= 0.) & (eta1_orig <= 1.) & (eta2_orig >= 0.) & (eta2_orig <= 1.))
    # Nothing to for particles that are in domain
    last_advec_percent[in_indices] = 100
    if (last_advec_percent == 100).all() :
        return


    where_out_e1_min = np.where(eta1_orig < 0.)
