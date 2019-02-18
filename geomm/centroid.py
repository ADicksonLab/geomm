import numpy as np

def centroid(coords, weights=None):
    """Return the centroid (AKA center of geometry or when weights
    correspond to masses the 'center of mass (COM)') of a set of
    coordinates.

    Parameters
    ----------

    coords : arraylike
        The coordinates to find the centroid of.

    weights : arraylike, optional
        Optional weights for the coordinates to compute a weighted
        centroid (center of mass).
       (Default = None)

    """

    # no weights
    if weights is None:
        centroid = np.average(coords, axis=0)
    # use the weights
    else:
        centroid = np.average(coords, axis=0, weights=weights)

    return centroid
