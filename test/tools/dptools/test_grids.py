#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
import numpy as np
import dptools.grids as grids
import nose.tools

FLOAT_TOLERANCE = 1E-12


def test_grid_gridcoord_to_cartesian_3d():
    grid = grids.Grid(
        origin=[ 1.0, 2.0, -3.0 ],
        basis=[[ 0.1, 0.0, 0.0 ], [ -0.4, 0.1, 0.0 ], [ 0.2, -0.3, 0.5 ]],
        ranges=[[ 0, 7 ], [ -5, 10 ], [  3, 13 ]]
    )
    coords = grid.gridcoord_to_cartesian([[ -1, 2, 9 ], [ -2, 4, 18 ]])
    true_coords = np.array([[ 1.9, -0.5, 1.5 ], [ 2.8, -3.0, 6.0 ]])
    diff = np.max(np.abs(coords - true_coords))
    assert diff < FLOAT_TOLERANCE


def test_grid_cartesian_to_gridcoord_3d():
    grid = grids.Grid(
        origin=[ 1.0, 2.0, -3.0 ],
        basis=[[ 0.1, 0.0, 0.0 ], [ -0.4, 0.1, 0.0 ], [ 0.2, -0.3, 0.5 ]],
        ranges=[[ 0, 7 ], [ -5, 10 ], [  3, 13 ]]
    )
    gridcoords = grid.cartesian_to_gridcoord([ 1.9, -0.5, 1.5 ])
    true_gridcoords = np.array([ -1, 2, 9 ])
    assert np.all(gridcoords == true_gridcoords)


def test_grid_get_corners_gridcoord():
    grid = grids.Grid(
        origin=[ 1.0, 2.0, -3.0 ],
        basis=[[ 0.1, 0.0, 0.0 ], [ -0.4, 0.15, 0.0 ],
                     [ 0.2, -0.3, 0.5 ]],
        ranges=[[ 0, 7 ], [ -5, 10 ], [  3, 13 ]]
    )
    corners = grid.get_corners(coordtype=grids.GRID_COORD)
    true_corners = np.array([[ 0, -5, 3 ], [ 0, -5, 13 ],
                             [ 0, 10, 3 ], [ 0, 10, 13 ],
                             [ 7, -5, 3 ], [ 7, -5, 13 ],
                             [ 7, 10, 3 ], [ 7, 10, 13 ]])
    assert np.all(corners == true_corners)


def test_grid_get_corners_cartesian():
    grid = grids.Grid(
        origin=[ 1.0, 2.0, -3.0 ],
        basis=[[ 0.1, 0.0, 0.0 ], [ -0.4, 0.15, 0.0 ],
                     [ 0.2, -0.3, 0.5 ]],
        ranges=[[ 0, 7 ], [ -5, 10 ], [  3, 13 ]]
    )
    corners = grid.get_corners(coordtype=grids.CARTESIAN_COORD)
    true_corners = np.array([[ 3.6, 0.35, -1.5 ], [ 5.6, -2.65, 3.5 ],
                             [ -2.4, 2.6, -1.5 ], [ -0.4, -0.4, 3.5 ],
                             [ 4.3, 0.35, -1.5 ], [ 6.3, -2.65, 3.5 ],
                             [ -1.7, 2.6, -1.5 ], [ 0.3, -0.4, 3.5 ]])
    diff = np.max(np.abs(corners - true_corners))
    assert diff < FLOAT_TOLERANCE


def test_grid_get_subgrid_ranges():
    grid1 = grids.Grid(origin=[ 1.0, 2.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 4 ], [ -2, 8 ]])
    grid2 = grids.Grid(origin=[ 4.0, 3.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -2, 2 ], [ -2, 2 ]])
    subranges = grid1.get_subgrid_ranges(grid2)
    true_subranges = np.array([[ 0, 4 ], [ -1, 3 ]])
    assert np.all(subranges == true_subranges)


@nose.tools.raises(ValueError)
def test_grid_get_subgrid_ranges_incompatible_gridvecs():
    grid1 = grids.Grid(origin=[ 1.0, 2.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 4 ], [ -2, 8 ]])
    grid2 = grids.Grid(origin=[ 4.1, 3.0 ],
                       basis=[[ 1.0, 0.1 ], [ 1.0, 1.0 ]],
                       ranges=[[ -2, 2 ], [ -2, 2 ]])
    subranges = grid1.get_subgrid_ranges(grid2)


@nose.tools.raises(ValueError)
def test_grid_get_subgrid_ranges_incompatible_origin():
    grid1 = grids.Grid(origin=[ 1.0, 2.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 4 ], [ -2, 8 ]])
    grid2 = grids.Grid(origin=[ 4.1, 3.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -2, 2 ], [ -2, 2 ]])
    subranges = grid1.get_subgrid_ranges(grid2)


@nose.tools.raises(ValueError)
def test_grid_get_subgrid_ranges_subgrid_too_big():
    grid1 = grids.Grid(origin=[ 1.0, 2.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 4 ], [ -2, 8 ]])
    grid2 = grids.Grid(origin=[ 4.0, 3.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -2, 3 ], [ -2, 2 ]])
    subranges = grid1.get_subgrid_ranges(grid2)


def test_grid_get_intersection_grid():
    grid1 = grids.Grid(origin=[ 8.0, 8.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 6 ], [ -3, 7 ]])
    grid2 = grids.Grid(origin=[ 16.0, 6.0 ],
                       basis=[[ 1.0, 0.0 ], [ 0.0, 1.0 ]],
                       ranges=[[ -3, 7 ], [ -4, 8 ]])
    intersec = grid1.get_intersection_grid(grid2)
    assert np.all(np.abs(intersec.origin - grid1.origin) < FLOAT_TOLERANCE)
    assert np.all(np.abs(intersec.basis - grid1.basis) < FLOAT_TOLERANCE)
    assert np.all(intersec.ranges == [[ -1, 6 ], [ -3, 6 ]])


def test_grid_get_gridpoints_grid():
    grid = grids.Grid(origin=[ 8.0, 8.0 ],
                      basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                      ranges=[[ -2, 0 ], [ -1, 1 ]])
    gridpoints = grid.get_gridpoints(coordtype=grids.GRID_COORD)
    true_gridpoints = np.array([[ -2, -1 ], [ -2, 0 ], [ -1, -1 ], [ -1, 0 ]])
    assert np.all(gridpoints == true_gridpoints)


def test_grid_get_gridpoints_cartesian():
    grid = grids.Grid(origin=[ 8.0, 8.0 ],
                      basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                      ranges=[[ -2, 0 ], [ -1, 1 ]])
    gridpoints = grid.get_gridpoints(coordtype=grids.CARTESIAN_COORD)
    true_gridpoints = np.array([[ 5.0, 7.0 ], [ 6.0, 8.0 ],
                                [ 6.0, 7.0 ], [ 7.0, 8.0 ]])
    diff = np.max(abs(gridpoints - true_gridpoints))
    assert diff < FLOAT_TOLERANCE


def test_griddata_get_subgrid_dataview():
    grid1 = grids.Grid(origin=[ 1.0, 2.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 4 ], [ -2, 8 ]])
    data1 = np.arange(7 * 10).reshape(( 7, 10 ))
    griddata = grids.GridData(grid1, data1)
    grid2 = grids.Grid(origin=[ 4.0, 3.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -2, 2 ], [ -2, 2 ]])
    data2 = griddata.get_subgrid_dataview(grid2)
    true_data2 = griddata.data[3:7, 1:5]
    diff = np.max(np.abs(data2 - true_data2))
    assert diff < FLOAT_TOLERANCE


def test_griddata_get_subgrid_dataview_modifying_data():
    grid1 = grids.Grid(origin=[ 1.0, 2.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 4 ], [ -2, 8 ]])
    data1 = np.arange(7 * 10).reshape(( 7, 10 ))
    griddata = grids.GridData(grid1, data1)
    grid2 = grids.Grid(origin=[ 4.0, 3.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -2, 2 ], [ -2, 2 ]])
    data2 = griddata.get_subgrid_dataview(grid2)
    data2 *= 2
    diff = np.max(np.abs(data2 - griddata.data[3:7, 1:5]))
    assert diff < FLOAT_TOLERANCE


def test_griddata_get_value_gridcoord():
    grid1 = grids.Grid(origin=[ 1.0, 2.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 4 ], [ -2, 8 ]])
    data1 = np.arange(7 * 10).reshape(( 7, 10 ))
    griddata = grids.GridData(grid1, data1)
    gridcoords = np.array([[ -3, -2 ], [ 3, 7 ]])
    data = griddata.get_value(gridcoords)
    true_data = np.array([ 0, 69 ])
    assert np.all(data == true_data)


def test_griddata_get_interpolated_value_gridcoord_coord_outside():
    grid1 = grids.Grid(origin=[ 1.0, 2.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 4 ], [ -2, 8 ]])
    data1 = np.arange(7 * 10).reshape(( 7, 10 ))
    griddata = grids.GridData(grid1, data1)
    gridcoords = np.array([[ -3, -2 ], [ 3, 8 ]])
    data = griddata.get_interpolated_value_gridcoord(gridcoords)
    true_data = np.array([ 0, 0 ])
    assert np.all(data == true_data)


def test_griddata_get_interpolated_value_cartesian():
    grid1 = grids.Grid(origin=[ 1.0, 2.0 ],
                       basis=[[ 1.0, 0.0 ], [ 1.0, 1.0 ]],
                       ranges=[[ -3, 4 ], [ -2, 8 ]])
    data1 = np.arange(7 * 10).reshape(( 7, 10 ))
    griddata = grids.GridData(grid1, data1)
    cartesians = grid1.gridcoord_to_cartesian([[ -3, -2 ], [ 3, 8 ]])
    # Perturb cartesians. Obtained value must be the same as we would have get
    # for the unperturbed ones.
    cartesians += 0.1
    data = griddata.get_interpolated_value(cartesians,
                                           coordtype=grids.CARTESIAN_COORD)
    true_data = np.array([ 0, 0 ])
    assert np.all(data == true_data)
