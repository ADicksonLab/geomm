import numpy as np

from geomm.centering import center

def group_pair_traj(traj, side_length_list, member_a_idxs, member_b_idxs):
    """Calls group_pair for each frame in traj"""
    return [group_pair(coords, unitcell_side_lengths, member_a_idxs, member_b_idxs)
            for coords,unitcell_side_lengths in zip(traj,side_length_list)]

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
    recentered_coords = center_complex(grouped_coords, unitcell_side_lengths, complex_idxs)

    return recentered_coords

def center_complex(coords, unitcell_side_lengths, complex_idxs):
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

    # center the complex
    centered_coords = center(coords, complex_centroid)

    return centered_coords


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


