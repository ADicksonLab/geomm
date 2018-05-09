import numpy as np

from geomm.theobald_qcp import theobald_qcp
from geomm.centering import center
from geomm.centroid import centroid

def superimpose(ref_coords, coords, centered=True, idxs=None, rot_mat=False, rmsd=False, weights=None):
    """Superimpose a set of coordinates to reference coordinates using the
    Theobald-QCP method.

    Inputs:

      ref_coords :: the template coordinates that `coords` will be aligned to.

      coords :: the coordinates which will be aligned to the template
      coordinates and have the rotation matrix applied to them and
      returned.

      centered (optional) :: Default `True`, if False will center
      the coordinates to be superimposed to the origine to apply the
      rotation matrix.

      idxs (optional) :: Default `None`. If given will superimpose the
      coordinates based only on the alignment on this subset of atoms
      to the reference.

      rot_mat (optional) :: Default `False`. If True will return the
      rotation matrix computed from the call to Theobald-QCP.

      rmsd (optional) :: Default `False`. If True will return the rmsd
      computed from the call to Theobald-QCP. NOTE: this may not be
      exactly the same as the computation of the RMSD by the
      definition (i.e. the explicit equation implemented here in
      `geomm.rmsd.calc_rmsd`)

      weights (optional) :: give weights to the coordinates for a
      weighted centroid ('center of mass') used in translations to the
      origin and to the centroid of the reference structure.

    """

    # first perform the theobald_qcp method to get the rotation matrix
    qcp_rmsd, rotation_matrix = theobald_qcp(ref_coords, coords,
                                             idxs=idxs,
                                             rot_mat=True, weights=weights)

    # if the coordinates are not already centered we must center them
    if not centered:
        centered_coords = center(coords, idxs=idxs, weights=weights)
    else:
        centered_coords = coords

    # rotate coords according to the rotation matrix
    rot_coords = np.dot(centered_coords, rotation_matrix)

    # find the centroid of the reference coordinates
    if idxs is not None:
        ref_centroid = centroid(ref_coords[idxs], weights=weights)
    else:
        ref_centroid = centroid(ref_coords, weights=weights)

    # translate the rotated coordinates to the reference centroid
    sup_coords = rot_coords + ref_centroid

    # return results based on what was asked for
    if rot_mat and rmsd:
        return rot_coords, rotation_matrix, qcp_rmsd
    elif rot_mat:
        return rot_coords, rotation_matrix
    elif rmsd:
        return rot_coords, qcp_rmsd
    else:
        return rot_coords
