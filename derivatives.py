from __future__ import division
import numpy as np
from globals_variables import *
from connectivity import get_face, connectivity
from geometry import geo


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

            if (face == 0) or (face == 2):
                deriv_func = np.zeros(NPTS2)
            if (face == 1) or (face == 3):
                deriv_func = np.zeros(NPTS1)

            # Computation:
            if face == 0:
                if (v0 == npat) and (f0 == face) :
                    deriv_func = finite_diff_internal_fwd(data[npat], f0, NPTS2)
                else :
                    deriv_func = cubic_spline_approximation(data[npat], face, \
                                                   data[v0], f0, NPTS2)
            elif face == 1:
                if (v1 == npat) and (f1 == face) :
                    deriv_func = finite_diff_internal_fwd(data[npat], f1, NPTS2)
                else :
                    deriv_func = cubic_spline_approximation(data[npat], face, \
                                                   data[v1], f1, NPTS2)
            elif face == 2:
                if (v2 == npat) and (f2 == face) :
                    # deriv_func = np.zeros(NPTS2)
                    deriv_func = - finite_diff_internal_fwd(data[npat], \
                                                            face, NPTS2)
                else :
                    deriv_func = - cubic_spline_approximation(data[npat], \
                                                    face, data[v2], f2, NPTS2)
            elif face == 3:
                if (v3 == npat) and (f3 == face) :
                    # deriv_func = np.zeros(NPTS1)
                    deriv_func = - finite_diff_internal_fwd(data[npat], face, \
                                                           NPTS2)
                else :
                    deriv_fun = - cubic_spline_approximation(data[npat], face, \
                                                   data[v3], f3, NPTS2)
            else:
                import sys, os
                sys.exit("Error in compute_derivatives_bound():" + \
                         " unrecognizable face ID.")

            derivatives_pat.append(deriv_func)
        derivatives.append(derivatives_pat)
    return derivatives


def finite_diff_internal_fwd(data, face, NPTS) :
    fip0 = get_face(data, face, indx=0)
    fip1 = get_face(data, face, indx=1)
    fip2 = get_face(data, face, indx=2)
    fip3 = get_face(data, face, indx=3)
    fip4 = get_face(data, face, indx=4)
    d2f  = (-25./12.*fip0 + 4.*fip1 - 3.*fip2 + 4./3.*fip3 - 0.25*fip4) \
           * (NPTS - 1)
    return d2f


def cubic_spline_approximation(data1, face1, data2, face2, NPTS) :
    """
    Computation of slopes using the cubic spline approximation.
    This function is inspired by the article:
    "Hermite Spline Interpolation on Patches for a Parallel
     Solving of the Vlasov-Poisson Equation"
    by Nicolas Crouseilles, Guillaume Latu, Eric Sonnendrucker

    data1 : array containing the data on patch 1 (P1)
    face1 : index of the face of P1 where P1 meets P2
    data2 : array containing the data on patch 2 (P2)
    face2 : index of the face of P2 where P2 meets P1
    NPTS  : number of points on the face1 and face2.
    """

    h = 1/NPTS
    [wp, wm] = compute_weights(h)

    res = 0.
    for j in range(10):
        fjp = get_face(data1, face1, indx=j+1)
        fjm = get_face(data2, face2, indx=j+1)
        res = res + wp[j] * fjp + wm[j] * fjm

    return res

def compute_weights(h):
    """
    Computation of weights to compute slopes for interpolation.
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
