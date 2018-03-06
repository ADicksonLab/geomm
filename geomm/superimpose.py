import numpy as np

from geomm.rmsd import theobald_qcp

def superimpose(ref_coords, coords, rot_mat=False, rmsd=False, weights=None):

    # first perform the theobald_qcp method to get the rotation matrix
    qcp_rmsd, rotation_matrix = theobald_qcp(ref_coords, coords, rot_mat=True, weights=weights)

    # rotate coords according to the rotation matrix
    rot_coords = np.dot(coords, rotation_matrix)

    if rot_mat and rmsd:
        return rot_coords, rot_mat, qcp_rmsd
    elif rot_mat:
        return rot_coords, rot_mat
    elif rmsd:
        return rot_coords, qcp_rmsd
    else:
        return rot_coords
