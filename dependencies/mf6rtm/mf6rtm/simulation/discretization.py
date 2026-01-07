"""
This module contains the Discretization functions for various ModFlow grid types.
"""

from mf6rtm.simulation.mf6api import Mf6API
import math


def total_cells_in_grid(modflow_api: Mf6API) -> int:
    return math.prod(grid_dimensions(modflow_api))


def grid_dimensions(modflow_api: Mf6API) -> tuple[int, ...]:
    grid_type = modflow_api.grid_type.upper()
    try:
        return __DISCRETIZATION_FUNCTIONS[grid_type](modflow_api)
    except KeyError or NotImplementedError:
        raise ValueError(f"Grid type '{grid_type}' is not yet supported.")


def __not_supported(*args, **kargs):
    raise NotImplementedError("This grid type is not supported.")


def __dis(api: Mf6API) -> tuple[int, int, int]:
    """
    Returns the total number of grid cells from the structured rectangular
    layered grid specified in the Discretization (DIS) Package.
    """
    simulation = api.sim
    discretization = simulation.get_model(simulation.model_names[0]).dis
    nlay = discretization.nlay.get_data()
    nrow = discretization.nrow.get_data()
    ncol = discretization.ncol.get_data()
    return (nlay, nrow, ncol)


def __disv(api: Mf6API) -> tuple[int, int]:
    """
    Returns the total number of grid cells from the unstructured layered grid
    specified in the Discretization by Vertices (DISV) Package.
    """
    simulation = api.sim
    discretization = simulation.get_model(simulation.model_names[0]).disv
    nlay = discretization.nlay.get_data()
    ncpl = discretization.ncpl.get_data()
    return (nlay, ncpl)


__DISCRETIZATION_FUNCTIONS = {
    "DIS": __dis,
    "DISV": __disv,
    "DISU": __not_supported,
    "UNDEFINED": __not_supported,
}
"""
Dictionary mapping Modflow 6 grid type strings to their corresponding
discretization functions.

Keys
----
DIS : str
    Structured rectangular layered grid type defined by the Discretization
    (DIS) Package.
DISV : str
    Unstructured layered grid type defined by the Discretization by Vertices
    (DISV) Package.
DISU : str
    Unstructured flexible grid type defined by the  Unstructured Discretization
    (DISU) Package. Not supported.
UNDEFINED : str
    Undefined grid type. Not supported.

Values
------
function
    Function to calculate total cells in the grid based on
    the grid type.
"""
