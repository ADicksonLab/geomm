import numpy as np

def group_pair(coords, unitcell_side_lengths, member_a_idxs, member_b_idxs):
    """For a pair of group of coordinates (e.g. atoms) this moves member_b
    coordinates to the image of the periodic unitcell that minimizes
    the difference between the centers of geometry between the two
    members (e.g. a protein and ligand).

    """


    # take the difference between the average coords of each
    # molecule.
    unitcell_half_lengths = unitcell_side_lengths * 0.5

    # initialize a new array for the coordinates
    grouped_coords = np.copy(coords)

    # get the coordinates for the ligand and member_b separately
    member_a_coords = coords[member_a_idxs, :]
    member_b_coords = coords[member_b_idxs, :]

    # calculate the centroids for each member
    member_a_centroid = member_a_coords.mean(axis=0)
    member_b_centroid = member_b_coords.mean(axis=0)

    # calculate the difference between them from the given coordinates
    centroid_dist = member_a_centroid - member_b_centroid

    # now we need to move one relative to the other by unitcell
    # lengths in order to find one that makes the difference between
    # the centroids smallest

    # When the centroid_dist are larger than the unitcell half length in
    # both the positive and negative direction

    # The positive direction

    # this gives us a tuple of two arrays. THe first array is the
    # indices of the frames that satisfied the boolean expression. The
    # second array is the index of the dimension that satisfied the
    # boolean expression, i.e. 0,1,2 for x,y,z. For example if we have
    # a frame where member_b is outside the box on two dimensions we
    # could have (array([0,0]), array([0,2])) for frame 1 on both x
    # and z dimensions. and only (array([0]), array([0])) for outside
    # on the first frame x dimension.
    pos_idxs = np.where(centroid_dist > unitcell_half_lengths)[0]

    # we simply pair elements from each of those arrays (pairwise) to
    # iterate over them. Essentially reshaping the two arrays
    for dim_idx in pos_idxs:
        # change the ligand coordinates to recenter them by adding the
        # cube side length in that direction
        grouped_coords[member_b_idxs, dim_idx] = (member_b_coords[:, dim_idx] +
                                                unitcell_side_lengths[dim_idx])

    # the negative direction is the same as the positive direction
    # except for the boolean expression in the where command and we
    # subtract the cube side length instead of adding it
    neg_idxs = np.where(centroid_dist < -unitcell_half_lengths)[0]
    for dim_idx in neg_idxs:
        grouped_coords[member_b_idxs, dim_idx] = (member_b_coords[:, dim_idx] -
                                                  unitcell_side_lengths[dim_idx])


    return grouped_coords


def apply_rectangular_pbcs(coords, unitcell_side_lengths, center_point=(0., 0., 0.,)):
    """Apply rectangular Periodic Boundary Conditions (PBCs) given the
    lengths of the unitcell and a center point positions of the box in
    the coordinate space. The default for the center point is (0,0,0)
    which is the case for OpenMM MD frames but not other MD systems.

    """

    # check to make sure everything looks okay
    assert len(coords.shape) == 2, \
        "coordinates should be rank 2 array, "\
        "this function operates on individual frames not trajectories."
    assert coords.shape[1] == 3, "coordinates are not of 3 dimensions"
    assert len(center_point) == 3, "center point is not of 3 dimensions"
    assert len(unitcell_side_lengths) == 3, "Unitcell side lengths are not of dimension 3"

    # cast the center point to an array
    center_point = np.array(center_point)

    # Calculate half box sizes
    unitcell_half_lengths = unitcell_side_lengths * 0.5

    # initialize the coordinates to be wrapped
    wrapped_coords = np.copy(coords)

    # find coords which are outside the box in the positive direction
    pos_idxs = np.where(coords > center_point + unitcell_half_lengths)

    # Groups the frame_idx, atom_idx, and dim_idx
    pos_idxs = list(zip(pos_idxs[0], pos_idxs[1]))

    # Restrict particle coordinates to the simulation box
    for atom_idx, dim_idx in pos_idxs:
        wrapped_coords[atom_idx, dim_idx] = (coords[atom_idx, dim_idx] -
                                                 unitcell_side_lengths[dim_idx])

    # Find where coords are less than  the negative half box sizes
    neg_idxs = np.where(coords < center_point - unitcell_half_lengths)

    # Groups the fram_idx, atom_idx and dim_idx where they are greater
    # than half box sizes
    neg_idxs = list(zip(neg_idxs[0], neg_idxs[1]))

    # Restrict particle coordinates to the simulation box
    for atom_idx, dim_idx in neg_idxs:
        wrapped_coords[atom_idx, dim_idx] = (coords[atom_idx, dim_idx] +
                                                 unitcell_side_lengths[dim_idx])

    return wrapped_coords

def recenter_pair(coords, unitcell_side_lengths, member_a_idxs, member_b_idxs):
    """
    This method moves group of ligand and member_b to the given new
    center. The geometric mean is used to find the center point of
    ligand and receptor group.

    """

    # group the pair
    grouped_coords = group_pair(coords, unitcell_side_lengths,
                                member_a_idxs, member_b_idxs)

    # then recenter them as a complex
    complex_idxs = (member_a_idxs, member_b_idxs)
    recentered_coords = recenter_complex(grouped_coords, unitcell_side_lengths, complex_idxs)

    return recentered_coords

def recenter_complex(coords, unitcell_side_lengths, complex_idxs):
    """Recenter a periodic unitcell around a complex (a list of lists of
    atoms idxs for each member of the complex), by first computing the
    centroids of each member then computing the centroid of the
    centroids, and using that as the new center of the box.

    """

    # compute the centroids of each member in the complex
    member_centroids = []
    for member_idxs in complex_idxs:
        centroid = coords[member_idxs].mean(axis=0)
        member_centroids.append(centroid)
    member_centroids = np.array(member_centroids)

    # compute the centroid of the centroids
    complex_centroid = member_centroids.mean(axis=0)

    # apply the periodic boundary conditions around the complex
    # centroid to recenter
    recentered_coords = apply_rectangular_pbcs(coords, unitcell_side_lengths,
                                               center_point=complex_centroid)

    return recentered_coords


def group_complex(coords, unitcell_side_lengths, complexes_idxs):
    """Given several groups of atom idxs (say molecules) put them in the
    image which minimizes the distance between the individual
    centroids to the collective centroid.

    """

    raise NotImplementedError

    unitcell_half_lengths = unitcell_side_lengths * 0.5

    grouped_coords = np.copy(coords)

    # get the coordinates of each complex
    complex_coords = []
    for complex_idxs in complexes_idxs:
        complex_coords.append(coords[complex_idxs, :])

    return grouped_coords

import numpy as np

def traj_box_vectors_to_lengths_angles(traj_box_vectors):
    """Convert box vectors for multiple 'frames' (a 'trajectory') to box lengths and angles."""

    traj_unitcell_lengths = []
    for basis in traj_box_vectors:
        traj_unitcell_lengths.append(np.array([np.linalg.norm(frame_v) for frame_v in basis]))

    traj_unitcell_lengths = np.array(traj_unitcell_lengths)

    traj_unitcell_angles = []
    for vs in traj_box_vectors:

        angles = np.array([np.degrees(
                            np.arccos(np.dot(vs[i], vs[j])/
                                      (np.linalg.norm(vs[i]) * np.linalg.norm(vs[j]))))
                           for i, j in [(0,1), (1,2), (2,0)]])

        traj_unitcell_angles.append(angles)

    traj_unitcell_angles = np.array(traj_unitcell_angles)

    return traj_unitcell_lengths, traj_unitcell_angles

def box_vectors_to_lengths_angles(box_vectors):
    """Convert box vectors for a single 'frame' to lengths and angles."""

    # calculate the lengths of the vectors through taking the norm of
    # them
    unitcell_lengths = []
    for basis in box_vectors:
        unitcell_lengths.append(np.linalg.norm(basis))
    unitcell_lengths = np.array(unitcell_lengths)

    # calculate the angles for the vectors
    unitcell_angles = np.array([np.degrees(
                        np.arccos(np.dot(box_vectors[i], box_vectors[j])/
                                  (np.linalg.norm(box_vectors[i]) * np.linalg.norm(box_vectors[j]))))
                       for i, j in [(0,1), (1,2), (2,0)]])

    return unitcell_lengths, unitcell_angles
