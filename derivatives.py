from __future__ import division
import numpy as np
from globals_variables import *
from connectivity import get_face, connectivity
from geometry import geo


def compute_weights2d(h):
    """
    Computation of weights to compute slopes for interpolation.
    This function is inspired by the article:
    "Hermite Spline Interpolation on Patches for a Parallel
     Solving of the Vlasov-Poisson Equation"
    by Nicolas Crouseilles, Guillaume Latu, Eric Sonnendrucker

    h : step of discretization"""

    weights1 = np.zeros((5, 3))

    weights1[1, 0] = - 0.822695035460993 / 2.0
    weights1[1, 1] = - 0.0141843971631206
    weights1[1, 2] = 0.00709219858156028

    weights1[2, 0] = 0.230496453900709 / 2.0
    weights1[2, 1] = 0.0212765957446808
    weights1[2, 2] = - 0.00443262411347518

    weights1[3, 0] = - 0.0567375886524823 / 2.0
    weights1[3, 1] = - 0.0141843971631206
    weights1[3, 2] = 0.

    weights1[4, 0] = 0.00797872340425532 / 2.0
    weights1[4, 1] = 0.00354609929078014
    weights1[4, 2] = 0.000443262411347518

    return weights1/h



def compute_slopes(tab, list_patchs, jac, tstep=0.):
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
                eta2_slopes[npat, 0] = deriv_face
            elif face == 1 :
                eta1_slopes[npat, 0] = deriv_face
            elif face == 2 :
                eta2_slopes[npat, 1] = deriv_face
            elif face == 3 :
                eta1_slopes[npat, 1] = deriv_face
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

            # Initialization
            deriv_func = np.zeros(NPTS2)

            # Computation:
            if face == 0:
                #Bottom-Left neighbors value and face:
                [[v_bl], [f_bl]] = connectivity(v1, (f1+1)%4)
                #Bottom-Right neighbors value and face:
                [[v_br], [f_br]] = connectivity(v3, (f3-1)%4)

                deriv_func = cubic_spline_approximation2d(data[npat], face, \
                                                        data[v0], f0, # bottom
                                                        data[v1], (f1+1)%4, # left
                                                        data[v_bl], f_bl, #BL
                                                        data[v3], (f3-1)%4, # right
                                                        data[v_br], f_br, # BR
                                                        NPTS2)
            elif face == 1:
                # Bottom-left:
                [[v_bl], [f_bl]] = connectivity(v0, (f0-1)%4)
                # Upper-left
                [[v_ul], [f_ul]] = connectivity(v2, (f2+1)%4)
                deriv_func = cubic_spline_approximation2d(data[npat], face, \
                                                        data[v1], f1, # left
                                                        data[v0], (f0-1)%4, # bottom
                                                        data[v_bl], f_bl, #BL
                                                        data[v2], (f2+1)%4, # upper
                                                        data[v_ul], f_ul, # UL
                                                        NPTS1)
            elif face == 2:
                # Upper-Left neighbors value and face:
                [[v_ul], [f_ul]] = connectivity(v1, (f1-1)%4)
                # Upper-Right neighbors value and face:
                [[v_ur], [f_ur]] = connectivity(v3, (f3+1)%4)
                deriv_func = - cubic_spline_approximation2d(data[npat], face, \
                                                          data[v2], f2,# upper
                                                          data[v1], (f1-1)%4, # left
                                                          data[v_ul], f_ul, #UL
                                                          data[v3], (f3+1)%4, # right
                                                          data[v_ur], f_ur, # UR
                                                          NPTS2)
            elif face == 3:
                # Bottom-right:
                [[v_br], [f_br]] = connectivity(v0, (f0+1)%4)
                # Upper-right
                [[v_ur], [f_ur]] = connectivity(v2, (f2-1)%4)
                deriv_func = - cubic_spline_approximation2d(data[npat], face, \
                                                          data[v3], f3, # right
                                                          data[v0], (f0+1)%4, # bottom
                                                          data[v_br], f_br, #BR
                                                          data[v2], (f2-1)%4, # upper
                                                          data[v_ur], f_ur, # UR
                                                          NPTS1)

            else:
                import sys, os
                sys.exit("Error in compute_derivatives_bound():" + \
                         " unrecognizable face ID.")

            derivatives_pat.append(deriv_func)
        derivatives.append(derivatives_pat)
    return derivatives


def cubic_spline_approximation2d(data_u, face_u,
                                 data_b, face_b,
                                 data_ul, face_ul,
                                 data_bl, face_bl,
                                 data_ur, face_ur,
                                 data_br, face_br,
                                 NPTS) :
    """
    Computation of 2D slopes using the cubic spline approximation.
    This function is inspired by the article:
    "Hermite Spline Interpolation on Patches for a Parallel
     Solving of the Vlasov-Poisson Equation"
    by Nicolas Crouseilles, Guillaume Latu, Eric Sonnendrucker
    """

    h = 1./NPTS
    wp = compute_weights2d(h)

    res = np.zeros(NPTS1)
    for l in range(3):
        for k in range(1,5):
            # initialization:
            fmkml = np.zeros(NPTS1)
            fmkpl = np.zeros(NPTS1)
            fpkml = np.zeros(NPTS1)
            fpkpl = np.zeros(NPTS1)
            if l == 0 :
                fmkml = get_face(data_b, face_b, indx=k)
                fmkpl = get_face(data_b, face_b, indx=k) # copy, w0 = w0/2
                fpkml = get_face(data_u, face_u, indx=k)
                fpkpl = get_face(data_u, face_u, indx=k) # ok, w0 = w0/2
            else :
                fmkml[l:] = get_face(data_b, face_b, indx=k)[:-l]
                fmkpl[:-l] = get_face(data_b, face_b, indx=k)[l:]
                fpkml[l:] = get_face(data_u, face_u, indx=k)[:-l]
                fpkpl[:-l] = get_face(data_u, face_u, indx=k)[l:]
                for m in range(l):
                    fmkml[m] = get_face(data_br, face_br, indx=k)[m]
                    fmkpl[-m-1] = get_face(data_bl, face_bl, indx=k)[-m-1]
                    fpkml[m] = get_face(data_ur, face_ur, indx=k)[m]
                    fpkpl[-m-1] = get_face(data_ul, face_ul, indx=k)[-m-1]
            # updating result
            res += wp[k,l] * (fmkml - fpkpl) + wp[k,l] * (fmkpl - fpkml)
    return res
