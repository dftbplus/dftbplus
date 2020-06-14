#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''I/O routines for grids.'''

import numpy as np
from . import grids


def scalarvtk(fname, griddata, varname='var'):
    """Save a scalar field from a GridData object as a structured 3D VTK ASCII
    file.

    Args:
        griddata (GridData): object containing grid and scalar field
        fname (filename or file handle) : destination file
        var (string): variable name
    """

    if hasattr(fname, 'write'):
        fh = fname
        closefh = False
    elif isinstance(fname, str):
        fh = open(fname, 'w')
        closefh = True
    else:
        raise ValueError('Cannot open fname as file')

    ndim = griddata.grid.dimension

    # Verify constrains on basis and data (scalar field, oriented along x,y,z)
    if len(griddata.data.shape) != 3:
        raise ValueError('GridData object must represent a 3D scalar vector')
    if (np.multiply(griddata.grid.basis, np.eye(3))
            - griddata.grid.basis > grids.FLOAT_TOLERANCE).any():
        raise ValueError('GridData Grid basis must be oriented along xyz')

    fh.write('# vtk DataFile Version 3.0 \n')
    fh.write('vtk output \n')
    fh.write('ASCII\n')
    fh.write('DATASET STRUCTURED_POINTS\n')
    fh.write('DIMENSIONS ')
    for dim in range(ndim):
        fh.write('{} '.format(griddata.grid.shape[dim]))
    fh.write('\n')
    fh.write('SPACING ')
    for dim in range(ndim):
        fh.write('{} '.format(griddata.grid.basis[dim, dim]))
    fh.write('\n')
    fh.write('ORIGIN ')
    for dim in range(ndim):
        fh.write('{} '.format(griddata.grid.origin[dim]))
    fh.write('\n')
    fh.write('POINT_DATA ')
    fh.write('{}\n'.format(griddata.data.size))
    fh.write('SCALARS {}%0A float\n'.format(varname))
    fh.write('LOOKUP_TABLE default\n')
    # Note: vtk format wants x varying faster in A(x,y,z) notation,
    # which correspond to column-major order (left index varying faster)
    for value in griddata.data.flatten(order='F'):
        fh.write('{}\n'.format(value))

    if closefh:
        fh.close()


def cube(fname, griddata, header='dptools cube file'):
    """Save a scalar field from a Grid Data oject as a Gaussian cube file.
        (without any atomistic structure)

    Args:
        griddata (GridData): object containing grid and scalar field
        fname (filename or file handle) : destination file
        header (string): cube file header
    """

    if hasattr(fname, 'write'):
        fh = fname
        closefh = False
    elif isinstance(fname, str):
        fh = open(fname, 'w')
        closefh = True
    else:
        raise ValueError('Cannot open fname as file')

    ndim = griddata.grid.dimension
    if ndim != 3:
        raise ValueError('cube output only valid for 3D systems')

    # Verify constrains on basis and data (scalar field, oriented along x,y,z)
    if len(griddata.data.shape) != 3:
        raise ValueError('GridData object must represent a 3D scalar vector')
    if (np.multiply(griddata.grid.basis, np.eye(3))
            - griddata.grid.basis > grids.FLOAT_TOLERANCE).any():
        raise ValueError('GridData Grid basis must be oriented along xyz')

    fh.write('#{}\n'.format(header))
    fh.write('#\n')
    fh.write('{} {} {} {} \n'.format(1, griddata.grid.origin[0],
                                     griddata.grid.origin[1],
                                     griddata.grid.origin[2]))
    fh.write('{} {} {} {}\n'.format(griddata.grid.shape[0],
                                    griddata.grid.basis[0, 0], 0.0, 0.0))
    fh.write('{} {} {} {}\n'.format(griddata.grid.shape[1],
                                    0.0, griddata.grid.basis[1, 1], 0.0))
    fh.write('{} {} {} {}\n'.format(griddata.grid.shape[2],
                                    0.0, 0.0, griddata.grid.basis[2, 2]))
    fh.write('{} {} {} {} {} \n'.format(1, 0.0, 0.0, 0.0, 0.0))
    #Note: cube file wants z varying indexes faster, i.e. row major order
    for value in griddata.data.flatten(order='C'):
        fh.write('{} '.format(value))

    if closefh:
        fh.close()
