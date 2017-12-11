#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Test for the grids module'''

import unittest
import numpy as np
import dptools.grids as grids

FLOAT_TOLERANCE = 1E-12


class GridTests(unittest.TestCase):
    '''Tests of the grid objects'''

    def test_gridcoord_to_cartesian_3d(self):
        '''Grid to cartesian coordinate in 3D'''
        grid = grids.Grid(
            origin=[1.0, 2.0, -3.0],
            basis=[[0.1, 0.0, 0.0], [-0.4, 0.1, 0.0], [0.2, -0.3, 0.5]],
            ranges=[[0, 7], [-5, 10], [3, 13]])
        coords = grid.gridcoord_to_cartesian([[-1, 2, 9], [-2, 4, 18]])
        true_coords = np.array([[1.9, -0.5, 1.5], [2.8, -3.0, 6.0]])
        diff = np.max(np.abs(coords - true_coords))
        self.assertLess(diff, FLOAT_TOLERANCE)


    def test_cartesian_to_gridcoord_3d(self):
        '''Cartesian to grid coordinate in 3D'''
        grid = grids.Grid(
            origin=[1.0, 2.0, -3.0],
            basis=[[0.1, 0.0, 0.0], [-0.4, 0.1, 0.0], [0.2, -0.3, 0.5]],
            ranges=[[0, 7], [-5, 10], [3, 13]]
        )
        gridcoords = grid.cartesian_to_gridcoord([1.9, -0.5, 1.5])
        true_gridcoords = np.array([-1.0, 2.0, 9.0])
        diff = np.max(np.abs(gridcoords - true_gridcoords))
        self.assertLess(diff, FLOAT_TOLERANCE)


    def test_get_corners_gridcoord(self):
        '''Getting the corners of in grid coordinates'''
        grid = grids.Grid(
            origin=[1.0, 2.0, -3.0],
            basis=[[0.1, 0.0, 0.0], [-0.4, 0.15, 0.0],
                   [0.2, -0.3, 0.5]],
            ranges=[[0, 7], [-5, 10], [3, 13]])
        corners = grid.get_corners(coordtype=grids.GRID_COORD)
        true_corners = np.array([[0, -5, 3], [0, -5, 13],
                                 [0, 10, 3], [0, 10, 13],
                                 [7, -5, 3], [7, -5, 13],
                                 [7, 10, 3], [7, 10, 13]])
        self.assertTrue(np.all(corners == true_corners))


    def test_get_corners_cartesian(self):
        '''Getting the corners in Cartesian coordinates'''
        grid = grids.Grid(
            origin=[1.0, 2.0, -3.0],
            basis=[[0.1, 0.0, 0.0], [-0.4, 0.15, 0.0],
                   [0.2, -0.3, 0.5]],
            ranges=[[0, 7], [-5, 10], [3, 13]])
        corners = grid.get_corners(coordtype=grids.CARTESIAN_COORD)
        true_corners = np.array([[3.6, 0.35, -1.5], [5.6, -2.65, 3.5],
                                 [-2.4, 2.6, -1.5], [-0.4, -0.4, 3.5],
                                 [4.3, 0.35, -1.5], [6.3, -2.65, 3.5],
                                 [-1.7, 2.6, -1.5], [0.3, -0.4, 3.5]])
        diff = np.max(np.abs(corners - true_corners))
        self.assertLess(diff, FLOAT_TOLERANCE)


    def test_get_subgrid_ranges(self):
        '''Getting subranges'''
        grid1 = grids.Grid(origin=[1.0, 2.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 4], [-2, 8]])
        grid2 = grids.Grid(origin=[4.0, 3.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-2, 2], [-2, 2]])
        subranges = grid1.get_subgrid_ranges(grid2)
        true_subranges = np.array([[0, 4], [-1, 3]])
        self.assertTrue(np.all(subranges == true_subranges))


    def test_get_subgrid_ranges_incompatible_gridvecs(self):
        '''Getting subgrid ranges for incompatible gridvecs'''
        grid1 = grids.Grid(origin=[1.0, 2.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 4], [-2, 8]])
        grid2 = grids.Grid(origin=[4.1, 3.0],
                           basis=[[1.0, 0.1], [1.0, 1.0]],
                           ranges=[[-2, 2], [-2, 2]])
        with self.assertRaises(ValueError):
            subranges = grid1.get_subgrid_ranges(grid2)


    def test_get_subgrid_ranges_incompatible_origin(self):
        grid1 = grids.Grid(origin=[1.0, 2.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 4], [-2, 8]])
        grid2 = grids.Grid(origin=[4.1, 3.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-2, 2], [-2, 2]])
        with self.assertRaises(ValueError):
            subranges = grid1.get_subgrid_ranges(grid2)


    def test_get_subgrid_ranges_subgrid_too_big(self):
        grid1 = grids.Grid(origin=[1.0, 2.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 4], [-2, 8]])
        grid2 = grids.Grid(origin=[4.0, 3.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-2, 3], [-2, 2]])
        with self.assertRaises(ValueError):
            subranges = grid1.get_subgrid_ranges(grid2)


    def test_get_intersection_grid(self):
        grid1 = grids.Grid(origin=[8.0, 8.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 6], [-3, 7]])
        grid2 = grids.Grid(origin=[16.0, 6.0],
                           basis=[[1.0, 0.0], [0.0, 1.0]],
                           ranges=[[-3, 7], [-4, 8]])
        intersec = grid1.get_intersection_grid(grid2)
        self.assertTrue(
            np.all(np.abs(intersec.origin - grid1.origin) < FLOAT_TOLERANCE))
        self.assertTrue(
            np.all(np.abs(intersec.basis - grid1.basis) < FLOAT_TOLERANCE))
        self.assertTrue(np.all(intersec.ranges == [[-1, 6], [-3, 6]]))


    def test_get_gridpoints_grid(self):
        grid = grids.Grid(origin=[8.0, 8.0],
                          basis=[[1.0, 0.0], [1.0, 1.0]],
                          ranges=[[-2, 0], [-1, 1]])
        gridpoints = grid.get_gridpoints(coordtype=grids.GRID_COORD)
        true_gridpoints = np.array(
            [[-2, -1], [-2, 0], [-1, -1], [-1, 0]])
        self.assertTrue(np.all(gridpoints == true_gridpoints))


    def test_get_gridpoints_cartesian(self):
        grid = grids.Grid(origin=[8.0, 8.0],
                          basis=[[1.0, 0.0], [1.0, 1.0]],
                          ranges=[[-2, 0], [-1, 1]])
        gridpoints = grid.get_gridpoints(coordtype=grids.CARTESIAN_COORD)
        true_gridpoints = np.array([[5.0, 7.0], [6.0, 8.0],
                                    [6.0, 7.0], [7.0, 8.0]])
        diff = np.max(abs(gridpoints - true_gridpoints))
        self.assertLess(diff, FLOAT_TOLERANCE)


class GridDataTests(unittest.TestCase):

    def test_get_subgrid_dataview(self):
        grid1 = grids.Grid(origin=[1.0, 2.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 4], [-2, 8]])
        data1 = np.arange(7 * 10).reshape((7, 10))
        griddata = grids.GridData(grid1, data1)
        grid2 = grids.Grid(origin=[4.0, 3.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-2, 2], [-2, 2]])
        data2 = griddata.get_subgrid_dataview(grid2)
        true_data2 = griddata.data[3:7, 1:5]
        diff = np.max(np.abs(data2 - true_data2))
        self.assertLess(diff, FLOAT_TOLERANCE)


    def test_get_subgrid_dataview_modifying_data(self):
        grid1 = grids.Grid(origin=[1.0, 2.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 4], [-2, 8]])
        data1 = np.arange(7 * 10).reshape((7, 10))
        griddata = grids.GridData(grid1, data1)
        grid2 = grids.Grid(origin=[4.0, 3.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-2, 2], [-2, 2]])
        data2 = griddata.get_subgrid_dataview(grid2)
        data2 *= 2
        diff = np.max(np.abs(data2 - griddata.data[3:7, 1:5]))
        self.assertLess(diff, FLOAT_TOLERANCE)


    def test_get_value_gridcoord(self):
        grid1 = grids.Grid(origin=[1.0, 2.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 4], [-2, 8]])
        data1 = np.arange(7 * 10).reshape((7, 10))
        griddata = grids.GridData(grid1, data1)
        gridcoords = np.array([[-3, -2], [3, 7]])
        data = griddata.get_value(gridcoords)
        true_data = np.array([0, 69])
        self.assertTrue(np.all(data == true_data))


    def test_get_interpolated_value_gc_coord_outside(self):
        grid1 = grids.Grid(origin=[1.0, 2.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 4], [-2, 8]])
        data1 = np.arange(7 * 10).reshape((7, 10))
        griddata = grids.GridData(grid1, data1)
        gridcoords = np.array([[-3, -2], [3, 8]])
        data = griddata.get_interpolated_value_gc(gridcoords)
        true_data = np.array([0, 0])
        self.assertTrue(np.all(data == true_data))


    def test_get_interpolated_value_cartesian(self):
        grid1 = grids.Grid(origin=[1.0, 2.0],
                           basis=[[1.0, 0.0], [1.0, 1.0]],
                           ranges=[[-3, 4], [-2, 8]])
        data1 = np.arange(7 * 10).reshape((7, 10))
        griddata = grids.GridData(grid1, data1)
        cartesians = grid1.gridcoord_to_cartesian([[-3, -2], [3, 8]])
        # Perturb cartesians. Obtained value must be the same as we would have
        # get for the unperturbed ones.
        cartesians += 0.1
        data = griddata.get_interpolated_value(cartesians,
                                               coordtype=grids.CARTESIAN_COORD)
        true_data = np.array([0, 0])
        self.assertTrue(np.all(data == true_data))


if __name__ == '__main__':
    unittest.main()
