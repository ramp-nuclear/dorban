"""
tests for the polybox class
"""
import numpy as np
from pytest import fixture

from dorban.geometry.boundary_conditions import Void
from dorban.geometry.cartesian import Cartesian


@fixture
def example_polybox():
    return Cartesian([np.array([1, 3]), np.array([2, 4]), np.array([5, 7])],
                     [Void()] * 6).to_polybox()



def test_refinment(example_polybox):
    geo = example_polybox
    split = [[2, 2, 2]] * geo.cells
    geo = geo.refine_mesh(split)
    lengths=[np.array([0.5,1,2.5])]*8+[np.array([1.5,1,2.5])]*8+[np.array([0.5,2,2.5])]*8+\
            [np.array([1.5,2,2.5])]*8+[np.array([0.5,1,3.5])]*8+[np.array([1.5,1,3.5])]*8+\
            [np.array([0.5,2,3.5])]*8+[np.array([1.5,2,3.5])]*8
    assert np.all(np.array(geo.lengths)==np.array(lengths))
