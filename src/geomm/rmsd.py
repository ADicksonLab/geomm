import numpy as np

def calc_rmsd(ref_coords, coords, idxs=None):
    """Calculate the RMSD between the reference and query coordinates.

    Parameters
    ----------

    ref_coords : arraylike
        The refence coordinates that will be aligned to.

    coords : arraylike
        The coordinates that will be rotated to match ref_coords.

    idxs : arraylike of int
        Indices of the atoms in the coords to actually compute the
        RMSD for.

    """

    # make sure the coords are the same size
    assert len(ref_coords.shape) == len(coords.shape), "rank of coords are not the same"
    assert ref_coords.shape[0] == coords.shape[0], \
        "Number of coordinates are not the same"

    # make sure the number of dimensions is 3
    assert (ref_coords.shape[1] == 3) and (coords.shape[1] == 3), \
        "Number of dimensions are not the same"

    if idxs is None:
        idxs = np.arange(ref_coords.shape[0])


    rmsd = np.sqrt(np.sum(np.square(coords[idxs, :] - ref_coords[idxs, :]),
                          axis=(0,1))/idxs.shape[0])

    return rmsd

