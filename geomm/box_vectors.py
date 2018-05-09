import numpy as np

def box_vectors_to_lengths_angles(box_vectors):
    """Convert box vectors for a single 'frame' to lengths and angles."""

    # calculate the lengths of the vectors through taking the norm of
    # them
    unitcell_lengths = []
    for basis in box_vectors:
        unitcell_lengths.append(np.linalg.norm(basis))
    unitcell_lengths = np.array(unitcell_lengths)

    # calculate the angles for the vectors
    unitcell_angles = np.array([np.degrees(
                        np.arccos(np.dot(box_vectors[i], box_vectors[j])/
                                  (np.linalg.norm(box_vectors[i]) * np.linalg.norm(box_vectors[j]))))
                       for i, j in [(0,1), (1,2), (2,0)]])

    return unitcell_lengths, unitcell_angles
