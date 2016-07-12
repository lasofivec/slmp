import numpy as np
from globals_variables import *

# *************************************
#
# This version is made to only work
# with the domain of a circle and
# an external crown made out 4 patches
# with the description as follows :
#
# *************************************

# Creation liaisons entre les patchs :
list_faces_duplicated = []
list_faces_duplicata = []
for dict in geo.connectivity:
    list_faces_duplicated.append(dict['original'])
    list_faces_duplicata.append(dict['clone'])
# TODO: how to treat exterior boundary conditions?


if domain == SQUARE_2p_domain:
    list_faces_duplicata.append([1,2])
    list_faces_duplicated.append([0,0])
    geo.external_faces.pop(0)
    geo.external_faces.pop(-2)
if domain == SQUARE_2p_C0_domain:
    list_faces_duplicata.append([1,3])
    list_faces_duplicated.append([0,1])
    geo.external_faces.pop(1)
    geo.external_faces.pop(-1)
if domain == SQUARE_4p_domain:
    list_faces_duplicata.append([1,3])
    list_faces_duplicated.append([0,1])
    list_faces_duplicata.append([2,1])
    list_faces_duplicated.append([3,3])
    geo.external_faces.pop(1)
    geo.external_faces.pop(2)
    geo.external_faces.pop(-1)
    geo.external_faces.pop(-3)


def connectivity(npat, face):
    """
    Returns the associated patch and face number to the couple (npat, face) passed in
    parameter

    Notes:
        Usual faces/edges description:
        ____2______
        |         |
        1    0    3
        |____0____|
    Args:
        npat: Id of the patch
        face: Id of the edge/face

    Returns:

    """

    if np.size(npat) == 1:
        if [npat, face] in geo.external_faces:
            return [npat, face]
        count0 = list_faces_duplicated.count([npat, face])
        if count0 == 1:
            return list_faces_duplicata[list_faces_duplicated.index([npat, face])]
        else:
            return list_faces_duplicated[list_faces_duplicata.index([npat, face])]
    else:
        list_pats = []#np.zeros_like(npat)
        list_faces = []#np.zeros_like(npat)
        for i in range(np.size(npat)):
            [pat_i, face_i] = connectivity(npat[i], face)
            list_pats.append(pat_i)
            list_faces.append(face_i)
        return [list_pats, list_faces]


def get_neighbours(npat):
    """
    Return a list of the 4 patches neighbouring at each face. If the face is a piece
    of the domain's boundary (ie. there is not another patch on that face) then
    it gives the number of the patch.
    Examples:
        The following geometry
         ____2_________2_____
        |         |          |
        1    0   3|1   1     3
        |____0____|____0_____|

        will give: get_neighbours(0) = (0, 0, 0, 1) (the only actual neighbour is at face 3, and it's patch#1)
                   get_neighbours(1) = (1, 0, 1, 1) (the only actual neighbour is at face 1, and it's patch#0)

    Args:
        npat: Number of patch from which we need the neighbours

    Returns:
        List of the neighbouring patches indices.
    """
    f0 = -1
    f1 = -1
    f2 = -1
    f3 = -1

    lid = []
    for i in range(0, 4):
        lid.append(connectivity(npat, i)[0])

    [f0, f1, f2, f3] = lid
    return f0, f1, f2, f3


def get_face(mat, face, indx=0):
    # ******************************************
    # Gives back a vector with the data
    # of the matrix "mat" at the face "face"
    # "indx" is the shiffting value from face
    # ******************************************
    if face == 0:
        return np.copy(mat[:, indx])
    elif face == 1:
        return np.copy(mat[indx, :])
    elif face == 2:
        return np.copy(mat[:, -1 - indx])
    elif face == 3:
        return np.copy(mat[-1 - indx, :])
    else:
        import sys
        sys.exit("Error in get_face(): face number should be between 0 and 4")


def mapping_norms(list_patchs, jac):
    # ******************************************
    # Function to compute the norm of the mapping
    # at the borders of each patch
    # ******************************************
    # Inspired by : V.D. Liseikin, Grid Generation Methods
    # DOI 10.1007/978-90-481-2912-6_2
    norms = []
    # computes the norm to every face
    for npat in list_patchs:
        norms12 = []
        for face in range(0, 4):
            # getting the mapping's derivatives at the border:
            d1F1 = get_face(jac[npat, 0], face)
            d2F1 = get_face(jac[npat, 1], face)
            d1F2 = get_face(jac[npat, 2], face)
            d2F2 = get_face(jac[npat, 3], face)
            # computing the norms of the direct and inverse mappings.
            sqrtg = d1F1 * d2F2 + d1F2 * d2F1
            sqrtg[np.where(sqrtg == 0.)] = epsilon
            norm1 = np.sqrt(d1F1 * d1F1 + d2F1 * d2F1)
            norm2 = np.sqrt(d1F2 * d1F2 + d2F2 * d2F2)
            norminv1 = np.sqrt(d2F2 ** 2 + d2F1 ** 2) / sqrtg
            norminv2 = np.sqrt(d1F2 ** 2 + d1F1 ** 2) / sqrtg
            # To avoid dividing by 0 :
            norm1[np.where(norm1 == 0.)] = 1.
            norm2[np.where(norm2 == 0.)] = 1.
            sqrtg[np.where(sqrtg == 0.)] = 1.
            # Computing the grads of the mapping
            # (vectors perp to tangent but not normalized)
            grad_F1 = [d1F1, d2F1]
            grad_F2 = [d1F2, d2F2]
            grad_invF1 = [d2F2, -d2F1]
            grad_invF2 = [-d1F2, d1F1]
            # normalization of norm vector :
            grad_F1 = grad_F1 / norm1
            grad_F2 = grad_F2 / norm2
            grad_invF1 = grad_invF1 / norminv1
            grad_invF2 = grad_invF2 / norminv2
            # storing
            norms12.append([grad_F1, grad_F2])
        norms.append(norms12)
    return norms


def finite_diff_internal_fwd(data, face, NPTS) :
    fip0 = get_face(data, face, indx=0)
    fip1 = get_face(data, face, indx=1)
    fip2 = get_face(data, face, indx=2)
    fip3 = get_face(data, face, indx=3)
    fip4 = get_face(data, face, indx=4)
    d2f  = (-25./12.*fip0 + 4.*fip1 - 3.*fip2 + 4./3.*fip3 - 0.25*fip4) * (NPTS - 1)
    return d2f


def finite_diff_internal_ctr(data1, face1, data2, face2, NPTS) :
    # Let's compute f_{i+1} :
    fip1 = get_face(data1, face1, indx=1)
    # Let's compute f_{i+2} :
    fip2 = get_face(data1, face1, indx=2)
    # Let's compute f_{i-1} :
    fim1 = get_face(data2, face2, indx=1)
    # Let's compute f_{i-2} :
    fim2 = get_face(data2, face2, indx=2)
    # Approximation of df1/d2 with 4 points :
    d2f = (NPTS - 1.) * (8.0 * (fip1 - fim1) - (fip2 - fim2)) / 12.0
    return d2f


def finite_diff_external_ctr(data, face, data_down, face_down, \
                             data_up, face_up, NPTS, sense, \
                             dw_isdir=False, up_isdir=False):
    fip1 = np.zeros(NPTS)
    fip2 = np.zeros(NPTS)
    fim1 = np.zeros(NPTS)
    fim2 = np.zeros(NPTS)
    # Determining sense:
    corner_indx = sense * -1
    # We want to compute derivative on eta1
    # Let's compute f_{i+1} :
    val_edg = get_face(data, face)
    fip1[:-1] = val_edg[1:]
    fip1[-1] = get_face(data_down, face_down, indx=1)[corner_indx]
    # Let's compute f_{i+2} :
    fip2[:-2] = val_edg[2:]
    fip2[-2] = get_face(data_down, face_down, indx=1)[corner_indx]
    fip2[-1] = get_face(data_down, face_down, indx=2)[corner_indx]
    # Let's compute f_{i-1} :
    fim1[1:] = val_edg[:-1]
    fim1[0] = get_face(data_up, face_up, indx=1)[corner_indx]
    # Let's compute f_{i-2} :
    fim2[2:] = val_edg[:-2]
    fim2[0] = get_face(data_up, face_up, indx=1)[corner_indx]
    fim2[1] = get_face(data_up, face_up, indx=2)[corner_indx]
    # Approximation of df1/d1 with 4 points :
    d1f = (NPTS - 1.) * (8.0 * (fip1 - fim1) - (fip2 - fim2)) / 12.0

    #Treating boundary conditions:
    if up_isdir:
        d1f[0] = (-25./12.*val_edg[0] + 4.*val_edg[1] - 3.*val_edg[2] \
                  + 4./3.*val_edg[3] - 0.25*val_edg[4]) * (NPTS - 1.)
        d1f[1] = (-25./12.*val_edg[1] + 4.*val_edg[2] - 3.*val_edg[3] \
                  + 4./3.*val_edg[4] - 0.25*val_edg[5]) * (NPTS - 1.)
    if dw_isdir:
        d1f[-1] = -(-25./12.*val_edg[-1] + 4.*val_edg[-2] - 3.*val_edg[-3] \
                  + 4./3.*val_edg[-4] - 0.25*val_edg[-5]) * (NPTS - 1.)
        d1f[-2] = -(-25./12.*val_edg[-2] + 4.*val_edg[-3] - 3.*val_edg[-4] \
                  + 4./3.*val_edg[-5] - 0.25*val_edg[-6]) * (NPTS - 1.)
    return d1f


def compute_derivatives_bound(data, list_patchs):
    # ******************************************
    # Fonction to compute partial derivatives
    # at borders of each patch
    # ******************************************
    # initializations :
    derivatives = []

    for npat in list_patchs:
        derivatives_pat = []
        # getting patchs neighbors:
        [v0, f0] = connectivity(npat, 0)
        [v1, f1] = connectivity(npat, 1)
        [v2, f2] = connectivity(npat, 2)
        [v3, f3] = connectivity(npat, 3)

        for face in range(4):
            up_isdir = False
            dw_isdir = False
            if (face == 0) or (face == 2):
                d1f = np.zeros(NPTS2)
                d2f = np.zeros(NPTS2)
            if (face == 1) or (face == 3):
                d1f = np.zeros(NPTS1)
                d2f = np.zeros(NPTS1)
            # Computation:
            if face == 0:
                if (f3 == 3) and (npat == v3):
                    up_isdir = True
                if (f1 == 1) and (npat == v1):
                    dw_isdir = True
                d1f = finite_diff_external_ctr(data[npat], face, data[v3], f3, \
                                               data[v1], f1, NPTS2, False, \
                                               dw_isdir, up_isdir)
                if (v0 == npat) and (f0 == face) :
                    d2f = finite_diff_internal_fwd(data[npat], f0, NPTS2)
                else :
                    d2f = finite_diff_internal_ctr(data[npat], face, \
                                                   data[v0], f0, NPTS2)
            if face == 1:
                if (v1 == npat) and (f1 == face) :
                    d1f = finite_diff_internal_fwd(data[npat], f1, NPTS2)
                else :
                    d1f = finite_diff_internal_ctr(data[npat], face, \
                                                   data[v1], f1, NPTS2)
                if (f2 == 2) and (npat == v2):
                    up_isdir = True
                if (f0 == 0) and (npat == v0):
                    dw_isdir = True
                d2f = finite_diff_external_ctr(data[npat], face, data[v2], f2, \
                                               data[v0], f0, NPTS2, True, \
                                               dw_isdir, up_isdir)
            if face == 2:
                if (f3 == 3) and (npat == v3):
                    up_isdir = True
                if (f1 == 1) and (npat == v1):
                    dw_isdir = True
                d1f = finite_diff_external_ctr(data[npat], face, data[v3], f3, \
                                               data[v1], f1, NPTS2, True, \
                                               dw_isdir, up_isdir)
                if (v2 == npat) and (f2 == face) :
                    d2f = - finite_diff_internal_fwd(data[npat], face, NPTS2)
                else :
                    d2f = - finite_diff_internal_ctr(data[npat], face, \
                                                   data[v2], f2, NPTS2)
            if face == 3:
                if (v3 == npat) and (f3 == face) :
                    d1f = -finite_diff_internal_fwd(data[npat], face, NPTS2)
                else :
                    d1f = -finite_diff_internal_ctr(data[npat], face, \
                                                   data[v3], f3, NPTS2)
                if (f2 == 2) and (npat == v2):
                    up_isdir = True
                if (f0 == 0) and (npat == v0):
                    dw_isdir = True
                d2f = finite_diff_external_ctr(data[npat], face, data[v2], f2,\
                                               data[v0], f0, NPTS2, False, \
                                               dw_isdir, up_isdir)
            derivatives_pat.append([d1f, d2f])
        derivatives.append(derivatives_pat)
    return derivatives


def compute_derivatives_bound2(data, list_patchs):
    """ Computes the derivatives at the boundaries of patches.
    This function is an alternative to compute_derivatives_bound(...)
    but this function uses numpy's gradient to compute the derivatives.
    Args:
      data (float array): contains the value of a function f.
      list_patchs (int) : number of patches in the domain.
    Return:
      derivatives (float array): array of size [npatchs, 4, 2, NPTS]
      contains for each patch, for each face, for each variable, the
      value of the gradient.
    """

    derivatives = []

    for npat in list_patchs:
        pat_grad = np.gradient(data[npat], \
                               *(1. / NPTS1, 1. / NPTS2), \
                               edge_order=1)
        faces_grad = []
        for face in range(4):
            d1f = get_face(pat_grad[0], face)
            d2f = get_face(pat_grad[1], face)
            faces_grad.append([d1f, d2f])
        derivatives.append(faces_grad)

    return derivatives


def compute_slopes(tab, list_patchs, jac):
    """ ******************************************
    # Function that computes the slopes (gradients)
    # It uses both the mapping norms and
    # the derivatives at the borders of each patch
    # ****************************************** """
    # initialization
    eta1_slopes = np.zeros((np.size(list_patchs), 2, NPTS2))
    eta2_slopes = np.zeros((np.size(list_patchs), 2, NPTS1))
    data = np.copy(tab.reshape((np.size(list_patchs), NPTS1, NPTS2)))

    # getting derivatives and norms.
    # ie. d(f)/d(eta_#) and n = (d(F)/d(eta_#), d(F)/d(eta_#))
    derivatives = compute_derivatives_bound(data, list_patchs)
    norms = mapping_norms(list_patchs, jac)

    for npat in list_patchs:

        deriv_npat = derivatives[npat]
        norms_npat = norms[npat]

        for face in range(4):

            deriv_face = deriv_npat[face]

            norms_face = norms_npat[face]

            if (which_deriv == 3):
                gradF1 = norms_face[0]
                gradF2 = norms_face[1]

                if (face == 0):
                    eta2_slopes[npat, 0] = deriv_face[0] * gradF2[0] + deriv_face[1] * gradF2[1]
                elif (face == 1):
                    eta1_slopes[npat, 0] = deriv_face[0] * gradF1[0] + deriv_face[1] * gradF1[1]
                elif (face == 2):
                    eta2_slopes[npat, 1] = deriv_face[0] * gradF2[0] + deriv_face[1] * gradF2[1]
                elif (face == 3):
                    eta1_slopes[npat, 1] = deriv_face[0] * gradF1[0] + deriv_face[1] * gradF1[1]
            else :
                import sys, os
                sys.exit("Error in compute_slopes(): not defined for this type of derivation method")

    return eta1_slopes, eta2_slopes


def transform_patch_coordinates(eta_tab, face_tab):
    """
    Transforms a point at distance eta in the boundary of index face
    to the coordinates in the standard patch.
    For example (0) in face 3 will give (1,0), whereas in face 0 it will be (0,0).
    IMPORTANT: It needs that all patches have the same orientation !
    Args:
        eta_tab: coordinate of point in the edge "face"
        face_tab: index of the edge of the patch where the point is.

    Returns:
        [eta1_prime, eta2_prime]: the coordinates in the standard patch
    """
    list_eta1 = []
    list_eta2 = []
    if np.size(eta_tab)==1:
        face = face_tab

        if face == 0:
            return [eta_tab, 0.0]
        elif face == 1:
            return[eta_tab, 0.0]
        elif face == 2:
            return [eta_tab, 1.0]
        elif face == 3:
            return [1.0, eta_tab]
        else:
            import sys
            sys.exit("ERROR in transform_patch_coordinates(): wrong face index")
            return
    for ind in range(np.size(eta_tab)):
        face = face_tab[ind]
        eta  = eta_tab[ind]
        if face == 0:
            list_eta1.append(eta)
            list_eta2.append(0.0)
        elif face == 1:
            list_eta1.append(0.0)
            list_eta2.append(eta)
        elif face == 2:
            list_eta1.append(eta)
            list_eta2.append(1.0)
        elif face == 3:
            list_eta1.append(1.0)
            list_eta2.append(eta)
        else:
            import sys
            sys.exit("ERROR in transform_patch_coordinates(): wrong face index")
    return [list_eta1, list_eta2]


def transform_advection(advec_coef1, advec_coef2, face, list_faces):
    # to transform A_i (advection in patch i: P_i) to A_j (resp. P_j)
    # we need to perform the following computation:
    # J_j * J_i * A_i = A_j
    # where J_i and J_j are the jacobian matrix of P_i and P_j
    if np.size(list_faces) == 1:
        list_faces = [list_faces]
    for ind in range(np.size(list_faces)):
        if np.abs(list_faces[ind] - face) == 1:
            temp = advec_coef1[ind] + 0.
            advec_coef1[ind] = advec_coef2[ind]
            advec_coef2[ind] = -temp
        elif list_faces[ind] == face:
            print "adv coef 1 =", advec_coef1
            print "adv coef 2 =", advec_coef2
            advec_coef1[ind] = -advec_coef1[ind]
            advec_coef2[ind] = -advec_coef2[ind]
    return [advec_coef1, advec_coef2]

def transform_advection2(advec_coef1, advec_coef2, where_from, where_to, eta1_from, eta2_from, eta1_to, eta2_to):
    # to transform A_i (advection in patch i: P_i) to A_j (resp. P_j)
    # we need to perform the following computation:
    # J_j * J_i * A_i = A_j
    # where J_i and J_j are the jacobian matrix of P_i and P_j

    size_prbm = np.size(where_from)

    if (size_prbm == 1) :
        u = eta1_from
        v = eta2_from
        npat = where_from
        D = geo[npat].evaluate_deriv(u, v, nderiv=1)
        d1F1 = D[1, :, :, 0][0]
        d2F1 = D[2, :, :, 0][0]
        d1F2 = D[1, :, :, 1][0]
        d2F2 = D[2, :, :, 1][0]

        u = eta1_to
        v = eta2_to
        npat = where_to
        D = geo[npat].evaluate_deriv(u, v, nderiv=1)
        d1G1 = D[1, :, :, 0][0]
        d2G1 = D[2, :, :, 0][0]
        d1G2 = D[1, :, :, 1][0]
        d2G2 = D[2, :, :, 1][0]

        a1 = advec_coef1
        a2 = advec_coef2

        new_advec1 = (d1G1*d1F1 + d2G1*d1F2)*a1 + (d1G1*d2F1 + d2G1*d2F2)*a2
        new_advec2 = (d1G2*d1F1 + d2G2*d1F2)*a1 + (d1G2*d2F1 + d2G2*d2F2)*a2

        return [new_advec1, new_advec2]
    
    new_advec1 = np.zeros_like(advec_coef1)
    new_advec2 = np.zeros_like(advec_coef2)

    for ind in range(size_prbm):
        u = eta1_from[ind]
        v = eta2_from[ind]
        npat = where_from[ind]
        D = geo[npat].evaluate_deriv(u, v, nderiv=1)
        d1F1 = D[1, :, :, 0]
        d2F1 = D[2, :, :, 0]
        d1F2 = D[1, :, :, 1]
        d2F2 = D[2, :, :, 1]

        u = eta1_to[ind]
        v = eta2_to[ind]
        npat = where_to[ind]
        D = geo[npat].evaluate_deriv(u, v, nderiv=1)
        d1G1 = D[1, :, :, 0]
        d2G1 = D[2, :, :, 0]
        d1G2 = D[1, :, :, 1]
        d2G2 = D[2, :, :, 1]

        # [d1F1, d2F1, d1F2, d2F2] = jacobian(geo, where_from[ind], eta1_from[ind], eta2_from[ind])
        # [d1G1, d2G1, d1G2, d2G2] = jacobian(geo, where_to[ind], eta1_to[ind], eta2_to[ind])

        a1 = advec_coef1[ind]
        a2 = advec_coef2[ind]
        
        new_advec1[ind] = (d1G1*d1F1 + d2G1*d1F2)*a1 + (d1G1*d2F1 + d2G1*d2F2)*a2
        new_advec2[ind] = (d1G2*d1F1 + d2G2*d1F2)*a1 + (d1G2*d2F1 + d2G2*d2F2)*a2
        
    return [new_advec1, new_advec2]

