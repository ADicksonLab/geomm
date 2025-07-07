import numpy as np
import pytest
from geomm.grouping import group_pair

def test_group_pair_no_shift():
    # member_b centroid is close to member_a, no shift needed
    coords = np.array([
        [1.0, 1.0, 1.0],  # member_a
        [2.0, 2.0, 2.0],  # member_b
    ])
    unitcell = np.array([10.0, 10.0, 10.0])
    member_a_idxs = [0]
    member_b_idxs = [1]
    result = group_pair(coords, unitcell, member_a_idxs, member_b_idxs)
    np.testing.assert_allclose(result, coords)

def test_group_pair_positive_shift():
    # member_b centroid is too low in x, should be shifted by +unitcell[0]
    coords = np.array([
        [1.0, 1.0, 1.0],   # member_a
        [-8.0, 1.0, 1.0],   # member_b
    ])
    unitcell = np.array([10.0, 10.0, 10.0])
    member_a_idxs = [0]
    member_b_idxs = [1]
    expected = np.array([
        [1.0, 1.0, 1.0],
        [2.0, 1.0, 1.0],  # shifted by +10 in x
    ])
    result = group_pair(coords, unitcell, member_a_idxs, member_b_idxs)
    np.testing.assert_allclose(result, expected)

def test_group_pair_negative_shift():
    # member_b centroid is to high in -y, should be shifted by -unitcell[1]
    coords = np.array([
        [1.0, 0.0, 1.0],   # member_a
        [1.0, 8.0, 1.0],   # member_b
    ])
    unitcell = np.array([10.0, 10.0, 10.0])
    member_a_idxs = [0]
    member_b_idxs = [1]
    expected = np.array([
        [1.0, 0.0, 1.0],
        [1.0, -2.0, 1.0],  # shifted by -10 in y
    ])
    result = group_pair(coords, unitcell, member_a_idxs, member_b_idxs)
    np.testing.assert_allclose(result, expected)

def test_group_pair_shift_multiple_dims():
    # member_b centroid is far in +x and -z, should be shifted in both
    coords = np.array([
        [1.0, 1.0, 8.0],   # member_a
        [8.0, 1.0, 1.0],   # member_b
    ])
    unitcell = np.array([10.0, 10.0, 10.0])
    member_a_idxs = [0]
    member_b_idxs = [1]
    expected = np.array([
        [1.0, 1.0, 8.0],
        [-2.0, 1.0, 11.0],  # shifted by -10 in x, +10 in z
    ])
    result = group_pair(coords, unitcell, member_a_idxs, member_b_idxs)
    np.testing.assert_allclose(result, expected)

def test_group_pair_multiple_atoms_per_group():
    # Test with multiple atoms per group
    coords = np.array([
        [1.0, 1.0, 1.0],   # member_a
        [1.0, 2.0, 1.0],   # member_a
        [8.0, 1.0, 1.0],   # member_b
        [8.0, 2.0, 1.0],   # member_b
    ])
    unitcell = np.array([10.0, 10.0, 10.0])
    member_a_idxs = [0, 1]
    member_b_idxs = [2, 3]
    expected = np.array([
        [1.0, 1.0, 1.0],
        [1.0, 2.0, 1.0],
        [-2.0, 1.0, 1.0],  # shifted by -10 in x
        [-2.0, 2.0, 1.0],  # shifted by -10 in x
    ])
    result = group_pair(coords, unitcell, member_a_idxs, member_b_idxs)
    np.testing.assert_allclose(result, expected)