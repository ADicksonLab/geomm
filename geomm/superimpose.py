import numpy as np

from geomm.theobald_qcp import theobald_qcp
from geomm.centering import center
from geomm.centroid import centroid

def superimpose_traj(ref_coords, traj,
                centered=False, idxs=None, weights=None,
                rot_mat=False, rmsd=False):
    """Superimpose a whole trajectory onto a set of ref_coords
    (see superimpose for more info)"""

    return [superimpose(ref_coords, coords, centered=centered, idxs=idxs, weights=weights, rot_mat=rot_mat, rmsd=rmsd) for coords in traj]

def superimpose(ref_coords, coords,
                centered=False, idxs=None, weights=None,
                rot_mat=False, rmsd=False):
    """Superimpose a set of coordinates to reference coordinates using the
    Theobald-QCP method.

    Inputs:

      ref_coords :: the template coordinates that `coords` will be aligned to.

      coords :: the coordinates which will be aligned to the template
      coordinates and have the rotation matrix applied to them and
      returned.

      centered (optional) :: Default `False`, if False will center
      the coordinates to be superimposed to the origin to apply the
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

    # find the centroid of the reference coordinates
    if idxs is not None:
        ref_centroid = centroid(ref_coords[idxs], weights=weights)
    else:
        ref_centroid = centroid(ref_coords, weights=weights)

    # if the coordinates are not already centered we must center them
    if not centered:
        centered_coords = center(coords, idxs=idxs, weights=weights)
        centered_ref_coords = center(ref_coords, idxs=idxs, weights=weights)
    else:
        centered_coords = coords
        centered_ref_coords = ref_coords

    # perform the theobald_qcp method to get the rotation matrix
    qcp_rmsd, rotation_matrix = theobald_qcp(centered_ref_coords, centered_coords,
                                             idxs=idxs,
                                             rot_mat=True, weights=weights)

    # rotate coords according to the rotation matrix
    rot_coords = np.dot(centered_coords, rotation_matrix)

    # translate the rotated coordinates to the reference centroid
    sup_coords = rot_coords + ref_centroid

    # return results based on what was asked for
    if rot_mat and rmsd:
        return sup_coords, rotation_matrix, qcp_rmsd
    elif rot_mat:
        return sup_coords, rotation_matrix
    elif rmsd:
        return sup_coords, qcp_rmsd
    else:
        return sup_coords

# the following method contains portions of the software mdtraj which
# is distributed under the following license
##############################################################################
# MDTraj: A Python Library for Loading, Saving, and Manipulating
#         Molecular Dynamics Trajectories.
# Copyright 2012-2014 Stanford University and the Authors
#
# Authors: Peter Eastman, Robert McGibbon
# Contributors: Kyle A. Beauchamp, Matthew Harrigan, Carlos Xavier Hernandez
#
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
#
# Portions of this code originate from the OpenMM molecular simulation
# toolkit, copyright (c) 2012 Stanford University and Peter Eastman. Those
# portions are distributed under the following terms:
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
##############################################################################

def alt_superimpose(ref_coords, coords):
    """Returns the translation and rotation mapping mobile onto target.
    Parameters
    ----------
    mobile : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a `single` frame, to be aligned onto target.
    target : ndarray, shape = (n_atoms, 3)
        xyz coordinates of a `single` frame
    Returns
    -------
    translation : ndarray, shape=(3,)
        Difference between the centroids of the two conformations
    rotation : ndarray, shape=(3,3)
        Rotation matrix to apply to mobile to carry out the transformation.
    """

    # ensure_type(mobile, 'float', 2, 'mobile', warn_on_cast=False, shape=(None, 3))
    # ensure_type(target, 'float', 2, 'target', warn_on_cast=False, shape=(target.shape[0], 3))

    mu1 = ref_coords.mean(0)
    mu2 = coords.mean(0)

    ref_coords = ref_coords - mu1
    coords = coords - mu2

    correlation_matrix = np.dot(np.transpose(ref_coords), coords)
    V, S, W_tr = np.linalg.svd(correlation_matrix)
    is_reflection = (np.linalg.det(V) * np.linalg.det(W_tr)) < 0.0
    if is_reflection:
        V[:, -1] = -V[:, -1]
    rotation = np.dot(V, W_tr)

    translation = mu2 - mu1.dot(rotation)

    return translation, rotation
