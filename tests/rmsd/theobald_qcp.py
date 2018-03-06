"""
Sample code to use the routine for fast RMSD & rotational matrix calculation.
For the example provided below, the minimum least-squares RMSD for the two
7-atom fragments should be 0.719106 A.

    And the corresponding 3x3 rotation matrix is:

    [[ 0.72216358 -0.52038257 -0.45572112]
     [ 0.69118937  0.51700833  0.50493528]
     [-0.0271479  -0.67963547  0.73304748]]


"""

import numpy as np
from geomm.rmsd import theobald_qcp, calc_rmsd

# Setup coordinates

frag_a = np.array([
    [-2.803, -15.373, 24.556],
    [0.893, -16.062, 25.147],
    [1.368, -12.371, 25.885],
    [-1.651, -12.153, 28.177],
    [-0.440, -15.218, 30.068],
    [2.551, -13.273, 31.372],
    [0.105, -11.330, 33.567],
])

frag_b = np.array([
    [-14.739, -18.673, 15.040],
    [-12.473, -15.810, 16.074],
    [-14.802, -13.307, 14.408],
    [-17.782, -14.852, 16.171],
    [-16.124, -14.617, 19.584],
    [-15.029, -11.037, 18.902],
    [-18.577, -10.001, 17.996]
])


N = frag_a.shape[0]
# Calculate center of geometry
comA = np.sum(frag_a, axis=0)/N
comB = np.sum(frag_b, axis=0)/N

# Center each fragment
frag_a = frag_a - comA
frag_b = frag_b - comB

# align to a subselection of the frames
idxs = np.array([0,1,2,3])
rmsd, rot_mat = theobald_qcp(frag_a, frag_b, idxs=idxs)

print("theobald alignment RMSD: {}".format(rmsd))

# apply the rotation matrix and calculate the rmsd of only the
# alignemnt coordinates to check
frag_b_rot = np.dot(frag_b, rot_mat)
rot_rmsd = calc_rmsd(frag_b_rot, frag_a, idxs=idxs)
print("Rotation alignment RMSD: {}".format(rot_rmsd))

# calculate the RMSD of the whole set of coordinates for the given
# alignment
frag_b_rot = np.dot(frag_b, rot_mat)
rot_rmsd = calc_rmsd(frag_b_rot, frag_a)
print("Total rotation alignment RMSD: {}".format(rot_rmsd))
