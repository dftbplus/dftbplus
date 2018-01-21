#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''General grid and grid data representations'''

import numpy as np
import numpy.linalg as la

# Floating point tolerance for equality
FLOAT_TOLERANCE = 1E-12


GRID_COORD = 1
CARTESIAN_COORD = 2


class Grid:
    """Represents parallelepiped grid with arbitrary number of dimensions.

    Args:
        origin: Cartesian coordinates of the grid origin. Shape: (ndim,)
        basis: Array of row vectors spanning the grid basis.
            Shape: (ndim, ndim).
        ranges: Ranges for the grid vector repetitions. Shape: (ndim, 2).

    Attributes:
        origin: Cartesian coordinates of the grid origin. Shape: (ndim,)
        basis: Array of row vectors spanning the grid basis.
            Shape: (ndim, ndim)
        ranges: Ranges for grid vector repetitions along the axis.
            Shape: (ndim, 2)
        upper_bounds: Upper range bound (exclusive) for every axis.
            Shape (ndim,)
        lower_bounds: Lower range bound (inclusive) for every axis.
            Shape (ndim,)
        shape: Nr. of grid points along each dimension. Shape: (ndim,).
        dimension: Dimension of the grid.
    """

    def __init__(self, origin, basis, ranges):
        self.origin = np.array(origin, dtype=float)
        self.basis = np.array(basis, dtype=float)
        self.ranges = np.array(ranges, dtype=int)
        self.lower_bounds = self.ranges[:, 0]
        self.upper_bounds = self.ranges[:, 1]
        self.shape = self.upper_bounds - self.lower_bounds
        self.shape.shape = (-1,)
        self._invbasis = la.inv(basis)
        self.dimension = len(self.origin)


    def set_origin(self, origin):
        """Sets a new origin for the grid.

        Args:
            origin: Cartesian coordinates of the new origin. Shape (ndim,).
        """
        self.origin = np.array(origin, dtype=float)


    def gridcoord_to_cartesian(self, gridcoords):
        """Returns cartesian coordinates for given grid coordinates.

        Args:
            gridcoords: Grid coordinates. Shape (ndim,) or (-1, ndim).

        Returns:
            Corresponding cartesian coordinates.
        """
        coords = np.dot(gridcoords, self.basis)
        coords += self.origin
        return coords


    def cartesian_to_gridcoord(self, cartcoords):
        """Returns grid coordinates for given Cartesian coordinates.

        Args:
            cartcoords: Cartesian coordinates. Shape (ndim,) or (-1, ndim)

        Returns:
            Corresponding grid coordinates.
        """
        gridcoords = np.dot(cartcoords - self.origin, self._invbasis)
        return gridcoords


    def get_corners(self, coordtype):
        """Returns the corners of the grid.

        Args:
            coordtype: Type of the returned coordinates, GRID_COORD and
                CARTESIAN_COORD for grid coordinates and Cartesian
                coordinates, respectively.

        Returns:
            Coordinates of the 2**ndim corners of the parallelepipedon.
        """
        # Gives tuples for upper lower range indices for each dimension
        # e.g. for 3D [[0, 0, 0], [0, 0, 1], [0, 1, 0], ..., [1, 1, 1]]
        corner_inds = np.indices((2,) * self.dimension)
        corner_inds = corner_inds.reshape((self.dimension, -1)).transpose()

        # Calc. corners by taking the lower/upper range value for each dimension
        # according corner_inds
        corner_gridcoords = []
        for inds in corner_inds:
            tmp = [np.take(self.ranges[ii], (inds[ii],))[0]
                   for ii in range(self.dimension)]
            corner_gridcoords.append(tmp)
        corner_gridcoords = np.vstack(corner_gridcoords)
        if coordtype == GRID_COORD:
            return corner_gridcoords
        else:
            return self.gridcoord_to_cartesian(corner_gridcoords)


    def get_gridpoints(self, coordtype):
        """Returns all grid points of the grid.

        Args:
            coordtype: Type of returned coordinates, GRID_COORD and
                CARTESIAN_COORD for grid coordinates and Cartesian
                coordinates, respectively.

        Returns:
            Coordinates of all grid points. Shape (-1, ngrid).
        """
        meshranges = [range(*mrange) for mrange in self.ranges]
        meshgrid = np.meshgrid(*meshranges)
        gridcoords = np.transpose(meshgrid).reshape(-1, self.dimension)
        if coordtype == GRID_COORD:
            return gridcoords
        else:
            return self.gridcoord_to_cartesian(gridcoords)


    def get_subgrid_ranges(self, subgrid):
        """Returns the grid coordinate ranges for a subgrid.

        Args:
            subgrid: Subgrid to look for.

        Returns:
            Ranges of grid coordinates which correspond to the given subgrid.

        Raises:
            ValueError: If subgrid not compatible with grid or not contained
                fully in it.
        """
        diff = np.abs(subgrid.basis - self.basis)
        if np.any(diff > FLOAT_TOLERANCE):
            raise ValueError("Incompatible gridvectors in subgrid")
        suborig_real = np.dot(subgrid.origin - self.origin, self._invbasis)
        suborig_int = np.rint(suborig_real).astype(int)
        diff = np.abs(suborig_real - suborig_int)
        if np.any(diff > FLOAT_TOLERANCE):
            raise ValueError("Incompatible grid origins")
        lower_bounds = subgrid.lower_bounds + suborig_int
        upper_bounds = subgrid.upper_bounds + suborig_int
        contained = np.all(lower_bounds >= self.lower_bounds)
        contained = contained and np.all(upper_bounds <= self.upper_bounds)
        if not contained:
            raise ValueError("Subgrid not fully contained in grid.")
        return np.transpose([lower_bounds, upper_bounds])


    def get_intersection_grid(self, other):
        """Returns a grid which represents the intersection with an other grid.

        Args:
            other: Other grid to intersect with.

        Returns:
            Smallest grid containing the entire intersection area. It has
            the same origin and basis as self.
        """
        assert self.dimension == other.dimension
        other_corners_cart = other.get_corners(coordtype=CARTESIAN_COORD)
        other_corners_grid = self.cartesian_to_gridcoord(other_corners_cart)
        other_lower_bounds = np.min(other_corners_grid, axis=0)
        other_upper_bounds = np.max(other_corners_grid, axis=0)
        lower_bounds = np.max([self.lower_bounds, other_lower_bounds], axis=0)
        upper_bounds = np.min([self.upper_bounds, other_upper_bounds], axis=0)
        lower_bounds = np.floor(lower_bounds).astype(int)
        upper_bounds = np.ceil(upper_bounds).astype(int)
        gridranges = np.vstack((lower_bounds, upper_bounds)).transpose()
        intersection_grid = Grid(self.origin.copy(), self.basis.copy(),
                                 gridranges)
        return intersection_grid


    def has_subgrid(self, subgrid):
        """Checks whether a grid is contained as subgrid.

        Args:
            subgrid: Subgrid to look for

        Returns:
            True if subgrid is a true subgrid of current grid, False otherwise.
        """
        try:
            _ = self.get_subgrid_ranges(subgrid)
        except ValueError:
            return False
        return True


    def has_gridcoord(self, gridcoord):
        """Cheks whether a given position is within the boundaries of the grid.

        Args:
            gridcoord: Grid coordinates. Shape: (-1, ndim) or (ndim,)

        Returns:
            True if the grid coordinate is within the grid region.
        """
        gridcoords = np.reshape(gridcoord, (-1, self.dimension))
        is_inside = np.logical_and(self.lower_bounds <= gridcoords,
                                   gridcoords < self.upper_bounds)
        is_inside = np.logical_and.reduce(is_inside, axis=1)
        return is_inside



class GridData:
    """Connects grid with volumetric data.

    Args:
        grid: Volumetric grid.
        data: Volumetric data.
    """

    def __init__(self, grid, data):
        self.grid = grid
        self.data = data


    def set_grid_origin(self, origin):
        """Sets origin of the grid.

        Args:
            origin: Cartesian coordinates of the origin.
        """
        self.grid.set_origin(origin)


    def get_subgrid_dataview(self, subgrid):
        """Returns a view to the data array for a given subgrid.

        Args:
            subgrid: Subgrid, which must be compatible and fully contained.

        Returns:
            View object to the corresponding part of the data array.
        """
        subgrid_ranges = self.grid.get_subgrid_ranges(subgrid)
        relative_ranges = subgrid_ranges - self.grid.lower_bounds[:, np.newaxis]
        sliceobj = [slice(*relrange) for relrange in relative_ranges]
        return self.data[sliceobj]


    def get_value(self, gridcoords):
        """Returns the value for a given (exact) grid position.

        Args:
            gridcoords: Grid coordinates. Shape: (-1, ndim) or (ndim,)

        Returns:
            Data value for given position.
        """
        datainds = gridcoords - self.grid.lower_bounds
        # Make a tuple out of it to force basic indexing
        datainds = tuple(datainds.transpose())
        return self.data[datainds]


    def get_interpolated_value(self, coords, coordtype):
        """Returns an interpolated value for an arbitrary position.

        Args:
            coords: Position coordinates. Shape: (-1, ndim) or (ndim,).  The
                interpolation is done via the get_interpolated_value_gc()
                method.
            coordtype: Coordinate type (GRID_COORD or CARTESIAN_COORD)

        Returns:
            Interpolated data.

        """
        if coordtype == GRID_COORD:
            gridcoords = coords
        else:
            gridcoords = self.grid.cartesian_to_gridcoord(coords)
        return self.get_interpolated_value_gc(gridcoords)


    def get_interpolated_value_gc(self, gridcoords):
        """Returns an interpolated value for an arbitrary grid position.

        In the current implementation this is the crudest possible
        interpolation, where grid positions are truncated and the appropriate
        data value is returned. For positions outside of the grid, zero will
        be returned.

        Args:
            gridcoords: Positions as grid coordinates. Shape (-1, ndim) or
                (ndim,).
        """
        values = np.zeros(len(gridcoords), dtype=self.data.dtype)
        gridcoords = np.floor(gridcoords).astype(int)
        is_inside = self.grid.has_gridcoord(gridcoords)
        gridcoords_inside = gridcoords[is_inside]
        data_inside = self.get_value(gridcoords_inside)
        np.place(values, is_inside, data_inside)
        return values
