import numpy as np

from dorban.system import Core
from dorban.geometry.boundary_conditions import Boundary
from dorban.utils import op_face


def debug_discontinuity_factors(core: Core, maximal_allowed_value: float = 5.0) -> list[tuple[int, int, int]]:
    """
    This is a tool to debug discontinuity factors. Usually discontinuity factors are between 0.25 and 4, this being the
    value of the discontinuity factors between fuel and water in a PWR type core. If the discontinuity values are outside
    this range it usually means that there is a bug in their definition.
    As discontinuity factors are the most complicated piece of information when defining a core it is the most common place
    to make mistakes. If it seems that a core is broken its is a good idea to run this tool.

    The function returns a list of the problematic cells. If the list is empty then all is good.
    The elements of the list are tuples of cell number, face number and neighbor number. The discontinuity factors between
    the cell and its given neighbor across the given face are wrong.

    If no discontinuity factors are defined an empty list will be returned.

    Parameters
    ----------
    core: Core
     The dorban core.
    maximal_allowed_value: float
     The maximal allowed value of a discontinuity factor. default is 5.

    Returns
    -------
    list[tuple[int, int, int]]
     A list of the problematic cells. If the list is empty then all is good.
     The elements of the list are tuples of cell number, face number and neighbor number. The discontinuity factors between
     the cell and its given neighbor across the given face are wrong.
    """
    if not 'df' in core.current_calc.__dict__:
        return []
    bad_cells = []
    for cell, neigh in enumerate(core.geometry.neighbors):
        for face, n in enumerate(neigh):
            if not isinstance(n, Boundary):
                m = np.max(core.current_calc.df[cell][face] / core.current_calc.df[n][op_face(face)])
                if m > maximal_allowed_value:
                    bad_cells.append(tuple([cell, face, n]))
    return bad_cells
