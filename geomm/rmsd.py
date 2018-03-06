import numpy as np

from geomm.pyqcprot import CalcRMSDRotationalMatrix

def calc_rmsd(ref_coords, coords, idxs=None):

    # make sure the coords are the same size
    assert len(ref_coords.shape) == len(coords.shape), "rank of coords are not the same"
    assert ref_coords.shape[0] == coords.shape[0], \
        "Number of coordinates are not the same"

    # make sure the number of dimensions is 3
    assert (ref_coords.shape[1] == 3) and (coords.shape[1] == 3), \
        "Number of dimensions are not the same"

    if idxs is None:
        idxs = np.arange(ref_coords.shape[0])


    rmsd = np.sqrt(np.sum(np.square(coords[idxs, :] - ref_coords[idxs, :]),
                          axis=(0,1))/idxs.shape[0])

    return rmsd

def theobald_qcp(ref_coords, coords, rot_mat=True, weights=None):
    """Wrapper around the pyqcprot implementation of the Theobald-QCP
    method for the calculation of RMSD and the RMSD minimizing rotation
    matrix. This function just gives a more pythonic API to the
    function

    Arguments:

    ref_coords :: the refence coordinates that will be aligned to

    coords :: the coordinates that will be rotated to match ref_coords

    rot_mat :: if True will return both the rmsd and the rotation matrix

    weights :: if your coordinates are weighted (e.g. mass) this is an array of those weights

    Returns:

    if rot_mat is True:

       rmsd :: (float) the rmsd of the two sets of coordinates

    if rot_mat is False:

       rmsd, rotation_matrix :: the rotation matrix that minimizes the RMSD

    """

    # make sure the coords are the same size
    assert ref_coords.shape[0] == coords.shape[0], \
        "Number of coordinates are not the same"

    # make sure the number of dimensions is 3
    assert (ref_coords.shape[1] == 3) and (coords.shape[1] == 3), \
        "Number of dimensions are not the same"

    # the number of coordinates (atoms)
    n_coords = ref_coords.shape[0]

    # make sure the weights if given are the right size
    if weights is not None:
        assert weights.shape[0] == n_coords, \
            "Number of weights given does not match the number of coordinates"

    # make an empty rotation matrix to fill from the cython function,
    # this is always an array of length 9
    rotation_matrix = np.zeros((9,), dtype=np.float64)

    # reshape the coordinates because this implementation requires 3xN
    # arrays and not Nx3 arrays
    ref_coords.reshape( (3, n_coords))

    # calculate the rmsd, and the rotation matrix which is modified in-place
    rmsd = CalcRMSDRotationalMatrix(ref_coords, coords,
                                    n_coords, rotation_matrix, weights)

    # if the rotation matrix was asked for return it, otherwise don't
    if rot_mat:
        # reshape the rotation matrix to be 2D
        return rmsd, rotation_matrix.reshape( (3, 3) )
    else:
        return rmsd
