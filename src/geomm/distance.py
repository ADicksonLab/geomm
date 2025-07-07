import numpy as np
import scipy.spatial.distance as dist
from scipy.spatial import KDTree

def distance():
    return None

def minimum_distance(coordsA, coordsB):
    """Calculate the minimum distance between members of coordsA and coordsB.
    Uses a fast binary search algorithm of order N*log(N).

    Parameters
    ----------

    coordsA : arraylike, shape (Natoms_A, 3)
        First set of coordinates.

    coordsB : arraylike, shape (Natoms_B, 3)
        Second set of coordinates.

    """

    # make sure the number of dimensions is 3
    assert (coordsA.shape[1] == 3) and (coordsB.shape[1] == 3), \
        "Minimum distance expecting arrays of shape (N, 3)"

    tree = KDTree(coordsA)
    return(tree.query(coordsB)[0].min())
