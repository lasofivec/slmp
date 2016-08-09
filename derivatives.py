from __future__ import division
import numpy as np
from globals_variables import *
from connectivity import get_face, connectivity
from geometry import geo

def compute_weights(h):
    """ Computation of weights to compute slopes for interpolation.
    This function is inspired by the article:
    "Hermite Spline Interpolation on Patches for a Parallel
     Solving of the Vlasov-Poisson Equation"
    by Nicolas Crouseilles, Guillaume Latu, Eric Sonnendrucker

    h : step of discretization"""

    # Initialization
    w_plus  = np.zeros(10)
    w_minus = np.zeros(10)
    
    # First we define some parameters:
    alpha = 1. - 2/14**2
    beta = alpha - 2/14**2/alpha

    w_plus[0]  = 1/beta * (39/49/h - 1/14/alpha * 3/49/h)
    w_plus[1]  = 1/beta * (-3/14/h -3/14**3/alpha/h)
    w_plus[2]  = 1/beta * (3/49/h - 39/14**2/49/alpha/h)
    w_plus[3]  = 1/beta * -3/14**2/h
    w_plus[4]  = 1/beta * 39 /14**2/49/alpha/h
    w_plus[5]  = 1/beta * (-3/14**3/alpha/h + 1/14**4/alpha/12/h)
    w_plus[6]  = 1/beta * ( 3/14**2/49/h - 8/14**4/alpha/12/h)
    w_plus[7]  = 1/beta * -3/14**4/alpha/h
    w_plus[8]  = 1/beta * 8/14**4/alpha/12/h
    w_plus[9]  = 1/beta * -1/14**4/alpha/12/h

    w_minus = -w_plus

    return [w_plus, w_minus]

def compute_slopes(tab, list_patchs, jac):
    """
    Function that computes the slopes (gradients) of the function at 
    the boundaries of the patch.

    tab: contains the data from which we need to compute the slopes
    list_patchs: list of the patchs' IDs.
    jac: jacobian matrix on each patch.
    """
    
    # initialization
    eta1_slopes = np.zeros((np.size(list_patchs), 2, NPTS2))
    eta2_slopes = np.zeros((np.size(list_patchs), 2, NPTS1))
    data = np.copy(tab.reshape((np.size(list_patchs), NPTS1, NPTS2)))

    derivs = compute_derivatives_bound(data, list_patchs)
    for npat in list_patchs:

        deriv_npat = derivs[npat]

        for face in range(4):

            deriv_face = deriv_npat[face]

            if face == 0:
                eta2_slopes[npat, 0] = deriv_face[1]
            elif face == 1 :
                eta1_slopes[npat, 0] = deriv_face[0]
            elif face == 2 :
                eta2_slopes[npat, 1] = deriv_face[1]
            elif face == 3 :
                eta1_slopes[npat, 1] = deriv_face[0]
            else :
                import sys, os
                sys.exit("Error in compute_slopes():" + \
                         " unrecognizable face ID.")

    return eta1_slopes, eta2_slopes


def compute_derivatives_bound(data, list_patchs):
    """
    Computes partial derivatives at borders of each patch

    tab: contains the data from which we need to compute the slopes
    list_patchs: list of the patchs' IDs.
    """
    # initializations :
    derivatives = []

    for npat in list_patchs:
        
        derivatives_pat = []
        
        # getting patch's connectivity on each face:
        [[v0], [f0]] = connectivity(npat, 0)
        [[v1], [f1]] = connectivity(npat, 1)
        [[v2], [f2]] = connectivity(npat, 2)
        [[v3], [f3]] = connectivity(npat, 3)

        for face in range(4):

            # Flags to know if there are dirichlet BC
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


# def compute_derivatives_bound2(data, list_patchs):
#     """ Computes the derivatives at the boundaries of patches.
#     This function is an alternative to compute_derivatives_bound(...)
#     but this function uses numpy's gradient to compute the derivatives.
#     Args:
#       data (float array): contains the value of a function f.
#       list_patchs (int) : number of patches in the domain.
#     Return:
#       derivatives (float array): array of size [npatchs, 4, 2, NPTS]
#       contains for each patch, for each face, for each variable, the
#       value of the gradient.
#     """

#     derivatives = []

#     for npat in list_patchs:
#         pat_grad = np.gradient(data[npat], \
#                                *(1. / NPTS1, 1. / NPTS2), \
#                                edge_order=1)
#         faces_grad = []
#         for face in range(4):
#             d1f = get_face(pat_grad[0], face)
#             d2f = get_face(pat_grad[1], face)
#             faces_grad.append([d1f, d2f])
#         derivatives.append(faces_grad)

#     return derivatives

# def mapping_norms(list_patchs, jac):
#     # ******************************************
#     # Function to compute the norm of the mapping
#     # at the borders of each patch
#     # ******************************************
#     # Inspired by : V.D. Liseikin, Grid Generation Methods
#     # DOI 10.1007/978-90-481-2912-6_2
#     norms = []
#     # computes the norm to every face
#     for npat in list_patchs:
#         norms12 = []
#         for face in range(0, 4):
#             # getting the mapping's derivatives at the border:
#             d1F1 = get_face(jac[npat, 0], face)
#             d2F1 = get_face(jac[npat, 1], face)
#             d1F2 = get_face(jac[npat, 2], face)
#             d2F2 = get_face(jac[npat, 3], face)
#             # computing the norms of the direct and inverse mappings.
#             sqrtg = d1F1 * d2F2 + d1F2 * d2F1
#             sqrtg[np.where(sqrtg == 0.)] = epsilon
#             norm1 = np.sqrt(d1F1 * d1F1 + d2F1 * d2F1)
#             norm2 = np.sqrt(d1F2 * d1F2 + d2F2 * d2F2)
#             norminv1 = np.sqrt(d2F2 ** 2 + d2F1 ** 2) / sqrtg
#             norminv2 = np.sqrt(d1F2 ** 2 + d1F1 ** 2) / sqrtg
#             # To avoid dividing by 0 :
#             norm1[np.where(norm1 == 0.)] = 1.
#             norm2[np.where(norm2 == 0.)] = 1.
#             sqrtg[np.where(sqrtg == 0.)] = 1.
#             # Computing the grads of the mapping
#             # (vectors perp to tangent but not normalized)
#             grad_F1 = [d1F1, d2F1]
#             grad_F2 = [d1F2, d2F2]
#             grad_invF1 = [d2F2, -d2F1]
#             grad_invF2 = [-d1F2, d1F1]
#             # normalization of norm vector :
#             grad_F1 = grad_F1 / norm1
#             grad_F2 = grad_F2 / norm2
#             grad_invF1 = grad_invF1 / norminv1
#             grad_invF2 = grad_invF2 / norminv2
#             # storing
#             norms12.append([grad_F1, grad_F2])
#         norms.append(norms12)
#     return norms


def finite_diff_internal_fwd(data, face, NPTS) :
    fip0 = get_face(data, face, indx=0)
    fip1 = get_face(data, face, indx=1)
    fip2 = get_face(data, face, indx=2)
    fip3 = get_face(data, face, indx=3)
    fip4 = get_face(data, face, indx=4)
    d2f  = (-25./12.*fip0 + 4.*fip1 - 3.*fip2 + 4./3.*fip3 - 0.25*fip4) * (NPTS - 1)
    return d2f


def finite_diff_internal_ctr(data1, face1, data2, face2, NPTS) :

    h = 1/NPTS
    [wp, wm] = compute_weights(h)

    res = 0.
    for j in range(10):
        fjp = get_face(data1, face1, indx=j+1)
        fjm = get_face(data2, face2, indx=j+1)
        res = res + wp[j] * fjp + wm[j] * fjm
    
    # # Let's compute f_{i+1} :
    # fip1 = get_face(data1, face1, indx=1)
    # # Let's compute f_{i+2} :
    # fip2 = get_face(data1, face1, indx=2)
    # # Let's compute f_{i-1} :
    # fim1 = get_face(data2, face2, indx=1)
    # # Let's compute f_{i-2} :
    # fim2 = get_face(data2, face2, indx=2)
    # # Approximation of df1/d2 with 4 points :
    # res = (NPTS - 1.) * (8.0 * (fip1 - fim1) - (fip2 - fim2)) / 12.0
    return res


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
    val_edg = get_face(data, face)
    # Let's compute f_{i+1} :
    fip1[:-1] = val_edg[1:]
    # TODO: should indx = 1, 2 or 0, 1 ???
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

