import numpy as np

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
    """This method applies periodic boundary conditions to the given
    coordinates"""

    # check to make sure everything looks okay
    assert len(coords.shape) == 2, \
        "coordinates should be rank 2 array, "
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
    pos_idxs = np.where(coords > unitcell_half_lengths + center_point)

    # Groups the frame_idx, atom_idx, and dim_idx
    pos_idxs = list(zip(pos_idxs[0], pos_idxs[1]))

    # Restrict particle coordinates to the simulation box
    for atom_idx, dim_idx in pos_idxs:
        wrapped_coords[atom_idx, dim_idx] = (coords[atom_idx, dim_idx] -
                                                 unitcell_side_lengths[dim_idx])

    # Find where coords are less than  the negative half box sizes
    neg_idxs = np.where(coords < -unitcell_half_lengths)

    # Groups the fram_idx, atom_idx and dim_idx where they are greater
    # than half box sizes
    neg_idxs = list(zip(neg_idxs[0], neg_idxs[1]))

    # Restrict particle coordinates to the simulation box
    for atom_idx, dim_idx in neg_idxs:
        wrapped_coords[atom_idx, dim_idx] = (coords[atom_idx, dim_idx] +
                                                 unitcell_side_lengths[dim_idx])

    return wrapped_coords

def recenter_pair(coords, unitcell_side_lengths, member_a_idxs, member_b_idxs,
                         new_box_center=np.array([0.0, 0.0, 0.0])):
    """
    This method moves group of ligand and member_b to the given new
    center. The geometric mean is used to find the center point of
    ligand and receptor group.

    """

    # Groups ligand and receptor together
    grouped_coords = group_pair(coords, unitcell_side_lengths,
                                member_a_idxs, member_b_idxs)

    # Combines the indices of ligand and receptor
    lig_receptor_idxs = np.concatenate((ligand_idxs, receptor_idxs))

    # Find coordinate of ligand_receptor group
    lig_receptor_coords = new_coords[:, lig_receptor_idxs, :]

    # calculates the geometric center of ligand_receptor
    geometric_center = lig_receptor_coords.mean(axis=1)

    # Find the translation vector coordinates
    position_vectors = new_box_center - geometric_center

    # Translate  system to the new_box_center using position_vectors
    for fram_idx, x in enumerate(new_coords):
        x += position_vectors[fram_idx]

    # Apply periodic boundary conditions to tranlated coordinates
    new_coords = apply_pbc(new_coords, unitcell_side_lengths)

    return new_coords
