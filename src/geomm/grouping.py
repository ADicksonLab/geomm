import numpy as np

from geomm.centering import center

def group_complex(coords, complexes_idxs):
    """TODO
    """

    raise NotImplementedError

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

    Parameters
    ----------

    coords : arraylike
        The coordinate array of the particles you will be
        transforming.

    unitcell_side_lengths : arraylike of shape (3)
        The lengths of the sides of a rectangular unitcell.

    member_a_idxs : arraylike of int of rank 1
        Collection of the indices that define that member of the pair.

    member_b_idxs : arraylike of int of rank 1
        Collection of the indices that define that member of the pair.

    Returns
    -------

    grouped_coords : arraylike
        Transformed coordinates.


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
        # change the ligand coordinates to center them by adding the
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

