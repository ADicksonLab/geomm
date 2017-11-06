import numpy as np

def recenter_receptor_ligand(positions, unitcell_side_lengths, ligand_idxs=None, receptor_idxs=None):
    """Sam's implementation of the recentering for a protein-ligand
    complex with cubic periodic boundaries."""

    unitcell_half_lengths = unitcell_side_lengths * 0.5

    # get the coordinates for the ligand and receptor separately
    lig_coords = positions[:, ligand_idxs, :]
    receptor_coords = positions[:, receptor_idxs, :]

    # take the difference between the average positions of each
    # molecule.

    # This takes the mean for each spatial dimension between all atoms
    # between the two sets of coordinates (molecules) and takes the
    # difference between them. The result is a 2D array with shape
    # (n_frames, 3)
    mean_diffs = receptor_coords.mean(axis=1) - lig_coords.mean(axis=1)

    # if the coordinates of the ligand are beyond the bounds
    # recenter. We have to test in two directions. When the mean_diffs
    # are larger than the unitcell half length in both the positive
    # and negative direction

    # The positive direction

    # this gives us a tuple of two arrays. THe first array is the
    # indices of the frames that satisfied the boolean expression. The
    # second array is the index of the dimension that satisfied the
    # boolean expression, i.e. 0,1,2 for x,y,z. FOr example if we have
    # a frame where the ligand is outside the box on two dimensions we
    # could have (array([0,0]), array([0,2])) for frame 1 on both x
    # and z dimensions. and only (array([0]), array([0])) for outside
    # on the first frame x dimension.
    pos_idxs = np.where(mean_diffs > unitcell_half_lengths)
    # we simply pair elements from each of those arrays (pairwise) to
    # iterate over them. Essentially reshaping the two arrays
    pos_idxs = list(zip(pos_idxs[0], pos_idxs[1]))
    for frame_idx, dim_idx in pos_idxs:
        # change the ligand coordinates to recenter them by adding the
        # cube side length in that direction
        lig_coords[frame_idx, :, dim_idx] = lig_coords[frame_idx, :, dim_idx] + \
                                                unitcell_side_lengths[frame_idx, dim_idx]

    # the negative direction is the same as the positive direction
    # except for the boolean expression in the where command and we
    # subtract the cube side length instead of adding it
    neg_idxs = np.where(mean_diffs < -unitcell_half_lengths)
    neg_idxs = list(zip(neg_idxs[0], neg_idxs[1]))
    for frame_idx, dim_idx in neg_idxs:
        lig_coords[frame_idx, :, dim_idx] = lig_coords[frame_idx, :, dim_idx] - \
                                                unitcell_side_lengths[frame_idx, dim_idx]

    # make a new array of the coords
    new_coords = positions.copy()
    # replace only the ligand coordinates
    new_coords[:, ligand_idxs, :] = lig_coords

    return new_coords
