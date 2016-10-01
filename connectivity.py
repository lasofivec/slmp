import numpy as np
from globals_variables import *
from geometry import geo

# Creation liaisons entre les patchs :
list_faces_duplicated = []
list_faces_duplicata = []
for dict in geo.connectivity:
    list_faces_duplicated.append(dict['original'])
    list_faces_duplicata.append(dict['clone'])
# TODO: how to treat exterior boundary conditions?

# For periodic BC we do this by hand:
if domain == SQUARE_2p_domain:
    list_faces_duplicata.append([1,2])
    list_faces_duplicated.append([0,0])
    geo.external_faces.pop(0)
    geo.external_faces.pop(-2)
    list_faces_duplicata.append([0,1])
    list_faces_duplicated.append([0,3])
    geo.external_faces.pop(0)
    geo.external_faces.pop(0)
    list_faces_duplicata.append([1,1])
    list_faces_duplicated.append([1,3])
    geo.external_faces.pop(0)
    geo.external_faces.pop(0)
if domain == SQUARE_2p_C0_domain:
    list_faces_duplicata.append([1,0])
    list_faces_duplicated.append([0,1])
    geo.external_faces.pop(1)
    geo.external_faces.pop(-3)
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
    Returns the associated patch and face number to the couple (npat, face)
    passed in parameter

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
            return [[npat], [face]]
        count0 = list_faces_duplicated.count([npat, face])
        if count0 == 1:
            ind = list_faces_duplicated.index([npat, face])
            return [[list_faces_duplicata[ind][0]],
                    [list_faces_duplicata[ind][1]]]
        else:
            ind = list_faces_duplicata.index([npat, face])
            return [[list_faces_duplicated[ind][0]],
                    [list_faces_duplicated[ind][1]]]
    else:
        list_pats = []
        list_faces = []
        for i in range(np.size(npat)):
            [pat_i, face_i] = connectivity(npat[i], face)
            list_pats.append(pat_i[0])
            list_faces.append(face_i[0])
        return [list_pats, list_faces]


def get_neighbours(npat):
    """
    Return a list of the 4 patches neighbouring at each face. If the face is a
    piece of the domain's boundary (ie. there is not another patch on that face)
    then it gives the number of the patch.
    Examples:
        The following geometry
         ____2_________2_____
        |         |          |
        1    0   3|1   1     3
        |____0____|____0_____|

        will give: get_neighbours(0) = (0, 0, 0, 1) (the only actual neighbour
                       is at face 3, and it's patch#1)
                   get_neighbours(1) = (1, 0, 1, 1) (the only actual neighbour
                       is at face 1, and it's patch#0)

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
    """
    Gives back a vector with the data of the matrix "mat"
    at the face "face" with "indx" shift from face.

    Args
    mat: matrix containing the value of the function at a given
         patch.
    face: index of the face where we wish to get the data from.
    indx: Shift value, if we want the values directly on the face,
          there is no need to give this argument. Else it gives the
          shift (or displacement) value.
    """
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


def transform_patch_coordinates(eta_tab, face_tab):
    """
    Transforms a point at distance eta in the boundary of index face
    to the coordinates in the standard patch.
    For example (0) in face 3 will give (1,0), whereas in face 0 it
    will be (0,0).
    IMPORTANT: It needs that all patches have the same orientation !
    Args:
        eta_tab: coordinate of point in the edge "face"
        face_tab: index of the edge of the patch where the point is.

    Returns:
        [eta1_prime, eta2_prime]: the coordinates in the standard patch
    """
    list_eta1 = []
    list_eta2 = []

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


# def transform_advection(advec_coef1, advec_coef2, face, list_faces):
#     # to transform A_i (advection in patch i: P_i) to A_j (resp. P_j)
#     # we need to perform the following computation:
#     # J_j * J_i * A_i = A_j
#     # where J_i and J_j are the jacobian matrix of P_i and P_j
#     if np.size(list_faces) == 1:
#         list_faces = [list_faces]
#     for ind in range(np.size(list_faces)):
#         if np.abs(list_faces[ind] - face) == 1:
#             temp = advec_coef1[ind] + 0.
#             advec_coef1[ind] = advec_coef2[ind]
#             advec_coef2[ind] = -temp
#         elif list_faces[ind] == face:
#             print "adv coef 1 =", advec_coef1
#             print "adv coef 2 =", advec_coef2
#             advec_coef1[ind] = -advec_coef1[ind]
#             advec_coef2[ind] = -advec_coef2[ind]
#     return [advec_coef1, advec_coef2]

def transform_advection(advec_coef1, advec_coef2,
                        where_from, where_to,
                        eta1_from, eta2_from,
                        eta1_to, eta2_to):
    # to transform A_i (advection in patch i: P_i) to A_j (resp. P_j)
    # we need to perform the following computation:
    # J_j * J_i * A_i = A_j
    # where J_i and J_j are the jacobian matrix of P_i and P_j

    size_prbm = np.size(where_from)

    if (size_prbm == 1) :
        a1_temp = advec_coef1
        a2_temp = advec_coef2

        u = eta1_from
        v = eta2_from
        npat = where_from
        D = geo[npat].evaluate_deriv(u, v, nderiv=1)
        d1F1 = D[1, :, :, 0][0]
        d2F1 = D[2, :, :, 0][0]
        d1F2 = D[1, :, :, 1][0]
        d2F2 = D[2, :, :, 1][0]

        a1 = d1F1 * a1_temp + d2F1 * a2_temp
        a2 = d1F2 * a1_temp + d2F2 * a2_temp

        u = eta1_to
        v = eta2_to
        npat = where_to
        D = geo[npat].evaluate_deriv(u, v, nderiv=1)
        d1G1 = D[1, :, :, 0][0]
        d2G1 = D[2, :, :, 0][0]
        d1G2 = D[1, :, :, 1][0]
        d2G2 = D[2, :, :, 1][0]
        sqrt_g = d1G1*d2G2 - d1G2*d2G1
        sqrt_g[np.where(sqrt_g == 0.)] = epsilon

        new_advec1 = (d2G2 * a1 - d2G1 * a2) / sqrt_g
        new_advec2 = (d1G1 * a2 - d1G2 * a1) / sqrt_g

        return [new_advec1, new_advec2]

    new_advec1 = np.zeros_like(advec_coef1)
    new_advec2 = np.zeros_like(advec_coef2)

    for ind in range(size_prbm):
        a1_temp = advec_coef1[ind]
        a2_temp = advec_coef2[ind]

        u = eta1_from[ind]
        v = eta2_from[ind]
        npat = where_from[ind]
        D = geo[npat].evaluate_deriv(u, v, nderiv=1)
        d1F1 = D[1, :, :, 0]
        d2F1 = D[2, :, :, 0]
        d1F2 = D[1, :, :, 1]
        d2F2 = D[2, :, :, 1]

        a1 = d1F1 * a1_temp + d2F1 * a2_temp
        a2 = d1F2 * a1_temp + d2F2 * a2_temp

        u = eta1_to[ind]
        v = eta2_to[ind]
        npat = where_to[ind]
        D = geo[npat].evaluate_deriv(u, v, nderiv=1)
        d1G1 = D[1, :, :, 0]
        d2G1 = D[2, :, :, 0]
        d1G2 = D[1, :, :, 1]
        d2G2 = D[2, :, :, 1]
        sqrt_g = d1G1*d2G2 - d1G2*d2G1
        sqrt_g[np.where(sqrt_g == 0.)] = epsilon

        new_advec1[ind] = (d2G2 * a1 - d2G1 * a2) / sqrt_g
        new_advec2[ind] = (d1G1 * a2 - d1G2 * a1) / sqrt_g

    return [new_advec1, new_advec2]
