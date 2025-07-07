import numpy as np
import pytest
from geomm.rmsd import calc_rmsd

def test_rmsd_identical_coords():
    ref = np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0]])
    coords = np.copy(ref)
    assert calc_rmsd(ref, coords) == 0.0

def test_rmsd_simple_translation():
    ref = np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0]])
    coords = ref + 1.0
    expected = np.sqrt(np.sum((coords - ref)**2) / ref.shape[0])
    assert np.isclose(calc_rmsd(ref, coords), expected)

def test_rmsd_with_idxs():
    ref = np.array([[0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0]])
    coords = np.array([[0.0, 0.0, 0.0],
                       [2.0, 0.0, 0.0],
                       [0.0, 2.0, 0.0]])
    idxs = np.array([1, 2])
    expected = np.sqrt(np.sum((coords[idxs] - ref[idxs])**2) / idxs.shape[0])
    assert np.isclose(calc_rmsd(ref, coords, idxs), expected)

def test_rmsd_shape_assertion():
    ref = np.zeros((3, 3))
    coords = np.zeros((4, 3))
    with pytest.raises(AssertionError):
        calc_rmsd(ref, coords)

def test_rmsd_dimension_assertion():
    ref = np.zeros((3, 2))
    coords = np.zeros((3, 2))
    with pytest.raises(AssertionError):
        calc_rmsd(ref, coords)

def test_rmsd_rank_assertion():
    ref = np.zeros((3, 3))
    coords = np.zeros((3, 3, 1))
    with pytest.raises(AssertionError):
        calc_rmsd(ref, coords)