#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

'''
Interface module for the communication between DFTB+ and
Python via the foreign function C-library ctypes.
'''


import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np
import numpy.linalg as la


# DFTB+ conversion factors
# (according to prog/dftb+/lib_common/constants.F90)
HARTREE__EV = 27.2113845
BOHR__AA = 0.529177249


class DftbPlus:
    '''ctypes interface to DFTB+.

    A shared library of DFTB+ is loaded and calculations of several
    physical quantities can be carried out. After completing the calculations,
    the results of the specified properties can be extracted easily.
    '''

    def __init__(self, libpath='./libdftbplus.so', hsdpath='./dftb_in.hsd',
                 logfile=None):
        '''Initializes a ctypes DFTB+ calculator object.

        Args:
            libpath (str): path to DFTB+ shared library
            hsdpath (str): path to DFTB+ input file
            logfile (str): name of log file
        '''

        self._libpath = libpath
        self._hsdpath = hsdpath
        self._logfile = logfile

        self._coords = None
        self._natoms = 0
        self._latvecs = None
        self._invlatvecs = None
        self._periodic = self._latvecs is not None
        self._relcoords = None
        self._props = None

        self._extpot = None
        self._extpotgrad = None

        # DFTB+ shared library
        self._dftbpluslib = ctypes.CDLL(self._libpath)

        input_str = ctypes.create_string_buffer(str.encode(self._hsdpath))

        # pointer to DFTB+ instance
        self._dftb_handler = ctypes.c_void_p()
        # pointer to DFTB+ input
        self._dftb_input = ctypes.c_void_p()

        # ndpointer instances used to describe 1darray
        # and 2darray in restypes and argtypes specifications
        self._dp1d = ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
        self._dp2d = ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')

        # Register supported (DFTB+) routines
        self._setup_interface()

        if self._logfile is not None:
            output_str = ctypes.create_string_buffer(str.encode(self._logfile))
            self._dftbpluslib.dftbp_init(self._dftb_handler, output_str)
        else:
            self._dftbpluslib.dftbp_init(self._dftb_handler, None)

        self._dftbpluslib.dftbp_get_input_from_file(self._dftb_handler,
                                                    input_str, self._dftb_input)

        self._dftbpluslib.dftbp_process_input(self._dftb_handler,
                                              self._dftb_input)

        self._natoms = self._dftbpluslib.dftbp_get_nr_atoms(self._dftb_handler)

        self._gradients = np.zeros((self._natoms, 3))
        self._energy = np.zeros(1)
        self._grosschg = np.zeros(self._natoms)


    def set_geometry(self, coords, latvecs=None, relcoords=False):
        '''Sets up the desired geometry.

        Args:
            coords (2darray):    absolute atomic positions
                                 (in atomic units)
            latvecs (2darray):   lattice vectors (in atomic units)
                                 (None for non-periodic structures)
            relcoords (2darray): relative atomic positions (in atomic units)
        '''

        self._coords = coords
        self._natoms = len(self._coords)
        self._latvecs = latvecs
        self._periodic = latvecs is not None

        if self._periodic:
            self._latvecs = np.array(latvecs, dtype=float)
            self._invlatvecs = la.inv(self._latvecs)
            if relcoords:
                self._relcoords = coords
                self._coords = np.dot(self._relcoords, self._latvecs)
            else:
                self._coords = coords
                self._relcoords = np.dot(self._coords, self._invlatvecs)
        else:
            self._latvecs = None
            self._invlatvecs = None
            self._relcoords = None

        if self._periodic:
            self._dftbpluslib.dftbp_set_coords_and_lattice_vecs(
                self._dftb_handler, self._coords, self._latvecs)
        else:
            self._dftbpluslib.dftbp_set_coords(self._dftb_handler, self._coords)


    def set_ext_pot(self, extpot, extpotgrad=None):
        '''Sets up an external potential.

        Args:
            extpot (2darray):     External potential at the position of each
                                  atom. Shape: [natom]. (in atomic units)
            extpotgrad (2darray): Gradient of the external potential at each
                                  atom. Shape: [natom, 3]. (in atomic units)
                                  This parameter is optional, you can pass None
                                  if you did not ask DFTB+ to calculate forces.
        '''

        self._extpot = extpot

        if extpotgrad is not None:
            self._extpotgrad = np.ctypeslib.as_ctypes(extpotgrad)

        self._dftbpluslib.dftbp_set_external_potential(
            self._dftb_handler, self._extpot, self._extpotgrad)


    def get_energy(self):
        '''Performs the energy calculation and queries
           the energy of the current geometry.

        Returns:
            energy (1darray): calculated energy (in atomic units)
        '''

        self._dftbpluslib.dftbp_get_energy(self._dftb_handler, self._energy)

        energy = self._energy[0]

        return energy


    def get_forces(self):
        '''Performs the calculation of the atomic forces and
           queries the forces of the current geometry.

        Returns:
            forces (2darray): calculated forces (in atomic units)
        '''

        self._dftbpluslib.dftbp_get_gradients(self._dftb_handler,
                                              self._gradients)

        forces = - self._gradients

        return forces


    def get_gross_charges(self):
        '''Queries the atomic Gross charges.

           Until state commit: 7c13358d (base: 19.1) of DFTB+:
           For non-trivial results, an energy or force calculation should be
           carried out beforehand. Otherwise, the atom populations are returned.
           More recent versions of DFTB+ adapt the behavior to the subroutines
           for the extraction of energy and gradients, therefore preventing the
           return of trivial charges without a calculation carried out.

           Sign convention: Electron has negative charge, so negative values
           indicate electron excess.

        Returns:
            grosschg (1darray): obtained Gross charges (in atomic units)
        '''

        self._dftbpluslib.dftbp_get_gross_charges(self._dftb_handler,
                                                  self._grosschg)

        grosschg = self._grosschg

        return grosschg


    def close(self):
        '''Finalizes the DFTB+ calculator.'''

        self._dftbpluslib.dftbp_final(self._dftb_handler)


    def _setup_interface(self):
        '''Bundle registrations of provided routines.
           Each time a wrap function is called to
           specify the argument and result types.
        '''

        self._wrap('dftbp_init', None,
                   [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p])

        self._wrap('dftbp_get_input_from_file', None,
                   [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p,
                    ctypes.POINTER(ctypes.c_void_p)])

        self._wrap('dftbp_process_input', None,
                   [ctypes.POINTER(ctypes.c_void_p),
                    ctypes.POINTER(ctypes.c_void_p)])

        self._wrap('dftbp_get_nr_atoms', ctypes.c_int,
                   [ctypes.POINTER(ctypes.c_void_p)])

        self._wrap('dftbp_set_coords', None,
                   [ctypes.POINTER(ctypes.c_void_p), self._dp2d])

        self._wrap('dftbp_set_coords_and_lattice_vecs',
                   None, [ctypes.POINTER(ctypes.c_void_p),
                          self._dp2d, self._dp2d])

        self._wrap('dftbp_set_external_potential', None,
                   [ctypes.POINTER(ctypes.c_void_p),
                    self._dp1d, ctypes.c_void_p])

        self._wrap('dftbp_get_energy', None,
                   [ctypes.POINTER(ctypes.c_void_p), self._dp1d])

        self._wrap('dftbp_get_gradients', None,
                   [ctypes.POINTER(ctypes.c_void_p), self._dp2d])

        self._wrap('dftbp_get_gross_charges', None,
                   [ctypes.POINTER(ctypes.c_void_p), self._dp1d])

        self._wrap('dftbp_final', None, [ctypes.POINTER(ctypes.c_void_p)])


    def _wrap(self, funcname, restype, argtypes):
        '''Wrap frequently used ctypes functions.

        Args:
            funcname (str):                     name of foreign function
            restype (C compatible data type):   ctypes type
                                                (result type of
                                                foreign function)
            argtypes (C compatible data types): tuple of ctypes types
                                                (argument types that
                                                foreign function accepts)
        '''

        func = self._dftbpluslib.__getattr__(funcname)
        func.restype = restype
        func.argtypes = argtypes
