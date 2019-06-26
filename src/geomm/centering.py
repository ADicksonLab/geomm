import numpy as np

from geomm.centroid import centroid

def center(coords, center_point):
    """Center coordinates at the origin based on a center point.

    If idxs are given the center of mass is computed only from those
    coordinates and if weights are given a weighted center of mass is
    computed.

    Parameters
    ----------

    coords : arraylike
       The coordinates you wish to center.

    center_point : arraylike of float of 1 dimension
       The point to center all other coordinates around. Should have
       one element for each dimension of coords

    Returns
    -------

    centered_coords : arraylike
        Transformed coordinates.

    """

    assert len(coords.shape) == 2, \
        "coordinates should be rank 2 array, "\
        "this function operates on individual frames not trajectories."
    assert coords.shape[1] == 3, "coordinates are not of 3 dimensions"
    assert len(center_point) == 3, "center point is not of 3 dimensions"

    return coords - center_point

def center_around(coords, idxs, weights=None):
    """Center coordinates at the origin based on a center point.

    If idxs are given the center of mass is computed only from those
    coordinates and if weights are given a weighted center of mass is
    computed.

    Parameters
    ----------

    coords : arraylike
       The coordinates you wish to center.

    idxs : arraylike of int
        The idxs of the coordinates to actually compute the centroid
        on, although the translation will act on all the coordinates.

    weights : arraylike of float, optional
        Give weights to the coordinates for a weighted centroid
        ('center of mass').

    Returns
    -------

    centered_coords : arraylike
        Transformed coordinates.

    """

    assert len(coords.shape) == 2, \
        "coordinates should be rank 2 array, "\
        "this function operates on individual frames not trajectories."
    assert coords.shape[1] == 3, "coordinates are not of 3 dimensions"
    assert len(idxs) > 0, "Must provide some idxs to compute a center of."

    return center(coords, centroid(coords[idxs], weights=weights))


def apply_rectangular_pbcs(coords, unitcell_side_lengths, center_point=(0., 0., 0.,)):
    """Apply rectangular Periodic Boundary Conditions (PBCs) given the
    lengths of the unitcell and a center point positions of the box in
    the coordinate space. The default for the center point is (0,0,0)
    which is the case for OpenMM MD frames but not other MD systems.

    Parameters
    ----------

    coords : arraylike
        The coordinate array of the particles you will be
        transforming.

    unitcell_side_lengths : arraylike of shape (3)
        The lengths of the sides of a rectangular unitcell.

    Returns
    -------

    wrapped_coords : arraylike
        Transformed coordinates. All fit within the box.

    Warning
    -------

    This method does not understand molecular topologies and will
    "break" bonds when moving molecules through boundaries.

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


def center_complex(coords, complex_idxs):
    """For a system with periodic boundary conditions move all members of
    a complex to the same image of the unitcell.

    Parameters
    ----------

    coords : arraylike
        The coordinate array of the particles you will be
        transforming.

    complex_idxs : list of arraylikes of int of rank 1
        A list where each member represents a member of the complex
        and is a collection of the indices that define that member.

    Returns
    -------

    centered_coords : arraylike
        Transformed coordinates.

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



