import numpy as np
import pytest
from geomm.distance import minimum_distance

def test_minimum_distance_simple():
    coordsA = np.array([[0, 0, 0], [1, 1, 1]])
    coordsB = np.array([[2, 2, 2], [3, 3, 3]])
    # Closest: [1,1,1] to [2,2,2] -> sqrt(3)
    expected = np.linalg.norm(np.array([1,1,1]) - np.array([2,2,2]))
    result = minimum_distance(coordsA, coordsB)
    assert np.isclose(result, expected)

def test_minimum_distance_identical_points():
    coordsA = np.array([[0, 0, 0], [1, 1, 1]])
    coordsB = np.array([[1, 1, 1], [2, 2, 2]])
    # [1,1,1] is in both, so min distance is 0
    assert minimum_distance(coordsA, coordsB) == 0.0

def test_minimum_distance_multiple_minima():
    coordsA = np.array([[0, 0, 0], [2, 2, 2]])
    coordsB = np.array([[1, 1, 1], [1, 1, 1]])
    # Both points in B are equidistant to both in A
    expected = np.linalg.norm(np.array([1,1,1]) - np.array([0,0,0]))
    result = minimum_distance(coordsA, coordsB)
    assert np.isclose(result, expected)

def test_minimum_distance_assertion_error():
    coordsA = np.array([[0, 0], [1, 1]])  # Not 3D
    coordsB = np.array([[2, 2], [3, 3]])
    with pytest.raises(AssertionError):
        minimum_distance(coordsA, coordsB)