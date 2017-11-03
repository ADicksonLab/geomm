import numpy as np

def recenter_receptor_ligand(positions, cube_side_lengths, ligand_idxs=None, receptor_idxs=None):
    """Sam's implementation of the recentering for a protein-ligand
    complex with cubic periodic boundaries."""

    unitcell_half_length = cube_side_lengths * 0.5

    # get the coordinates for the ligand and receptor separately
    lig_coords = positions[:, ligand_idxs, :]
    receptor_coords = positions[:, receptor_idxs, :]

    # take the difference between the average positions of each
    # molecule
    mean_diffs = receptor_coords.mean(axis=1) - lig_coords.mean(axis=1)

    # if the coordinates of the ligand are beyond the bounds recenter
    pos_idxs = np.where(mean_diffs > unitcell_half_length)
    pos_idxs = list(zip(pos_idxs[0], pos_idxs[1]))
    for pos_idx in pos_idxs:
        lig_coords[pos_idx[0], :, pos_idx[1]] = lig_coords[pos_idx[0], :, pos_idx[1]] + \
                                                cube_side_lengths[pos_idx]
    neg_idxs = np.where(mean_diffs < -unitcell_half_length)
    neg_idxs = list(zip(neg_idxs[0], neg_idxs[1]))
    for neg_idx in neg_idxs:
        lig_coords[neg_idx[0], :, neg_idx[1]] = lig_coords[neg_idx[0], :, neg_idx[1]] - \
                                                cube_side_lengths[neg_idx]

    # make a new array of the coords
    new_coords = positions.copy()
    new_coords[:, ligand_idxs, :] = lig_coords

    return new_coords
