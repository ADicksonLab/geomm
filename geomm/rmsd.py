import numpy as np

from geomm.pyqcprot import CalcRMSDRotationalMatrix

def rmsd(traj, ref, idx=None):
    # calculates rmsd between frames of two trajectories: traj and ref
    if len(idx) == 0:
        idx = range(len(traj))
    # check if traj and ref have the same number of atoms
    if len(traj[idx]) != len(ref[idx]):
        raise Exception('RMSD calculation must use the same number of atoms in each array')

    return np.sqrt(np.sum(np.square(traj[idx, :] - ref[idx, :]),
                          axis=(0,1))/idx.shape[0])


def theobald_qcp(coords_A, coords_B, rot_mat=True, weights=None):
    """Wrapper around the pyqcprot implementation of the Theobald-QCP
    method for the calculation of RMSD and the RMSD minimizing rotation
    matrix. This function just gives a more pythonic API to the
    function

    Arguments:

    coords_A :: the refence coordinates that will be aligned to

    coords_B :: the coordinates that will be rotated to match coords_A

    rot_mat :: if True will return both the rmsd and the rotation matrix

    weights :: if your coordinates are weighted (e.g. mass) this is an array of those weights

    Returns:

    if rot_mat is True:

       rmsd :: (float) the rmsd of the two sets of coordinates

    if rot_mat is False:

       rmsd, rotation_matrix :: the rotation matrix that minimizes the RMSD

    """

    # make sure the coords are the same size
    assert coords_A.shape[0] == coords_B.shape[0], "Number of coordinates are not the same"
    # make sure the number of dimensions is 3
    assert (coords_A.shape[1] == 3) and (coords_B.shape[1] == 3),\
        "Number of dimensions are not the same"

    # the number of coordinates (atoms)
    n_coords = coords_A.shape[0]

    # make sure the weights if given are the right size
    if weights is not None:
        assert weights.shape[0] == n_coords, \
            "Number of weights given does not match the number of coordinates"

    # make an empty rotation matrix to fill from the cython function,
    # this is always an array of length 9
    rotation_matrix = np.zeros((9,), dtype=np.float64)

    # reshape the coordinates because this implementation requires 3xN
    # arrays and not Nx3 arrays
    coords_A.reshape( (3, n_coords))

    # calculate the rmsd, and the rotation matrix which is modified in-place
    rmsd = CalcRMSDRotationalMatrix(coords_A, coords_B,
                                    n_coords, rotation_matrix, weights)

    # if the rotation matrix was asked for return it, otherwise don't
    if rot_mat:
        # reshape the rotation matrix to be 2D
        return rmsd, rotation_matrix.reshape( (3, 3) )
    else:
        return rmsd
