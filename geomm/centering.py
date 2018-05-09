import numpy as np

from geomm.centroid import centroid

def center(coords, idxs=None, weights=None):
    """Center coordinates at the origing based on their center of mass. If
    idxs are given the center of mass is computed only from those
    coordinates and if weights are given a weighted center of mass is
    computed.

    Inputs:

        coords :: the coordinates you wish to center

        idxs (optional) :: the idxs of the coordinates to actually
        compute the centroid on, although the translation will act on
        all the coordinates.

        weights (optional) :: give weights to the coordinates for a
           weighted centroid ('center of mass')

    """

    if idxs is None:
        centered_coords = coords - centroid(coords, weights=weights)
    else:
        centered_coords = coords - centroid(coords[idxs], weights=weights)

    return centered_coords

