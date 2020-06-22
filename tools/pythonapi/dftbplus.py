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


import os
import ctypes
from numpy.ctypeslib import ndpointer
import numpy as np


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


    def __init__(self, libpath='./libdftbplus', hsdpath='./dftb_in.hsd',
                 logfile=None):
        '''Initializes a ctypes DFTB+ calculator object.

        Args:
            libpath (str): path to DFTB+ shared library
            hsdpath (str): path to DFTB+ input file
            logfile (str): name of log file
        '''

        self._natoms = 0

        # DFTB+ shared library
        # use convenient Numpy wrapper to take different possible extensions
        # of the shared library into account (operating system dependent).
        libht = os.path.split(libpath)
        self._dftbpluslib = np.ctypeslib.load_library(libht[1], libht[0])

        input_str = ctypes.create_string_buffer(str.encode(hsdpath))

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

        if logfile is not None:
            output_str = ctypes.create_string_buffer(str.encode(logfile))
            self._dftbpluslib.dftbp_init(self._dftb_handler, output_str)
        else:
            self._dftbpluslib.dftbp_init(self._dftb_handler, None)

        self._dftbpluslib.dftbp_get_input_from_file(self._dftb_handler,
                                                    input_str, self._dftb_input)

        self._dftbpluslib.dftbp_process_input(self._dftb_handler,
                                              self._dftb_input)

        self._natoms = self._dftbpluslib.dftbp_get_nr_atoms(self._dftb_handler)


    def set_geometry(self, coords, latvecs=None):
        '''Sets up the desired geometry.

        Args:
            coords (2darray):    absolute atomic positions
                                 (in atomic units)
            latvecs (2darray):   lattice vectors (in atomic units)
                                 (None for non-periodic structures)
        '''

        periodic = latvecs is not None

        if periodic:
            latvecs = np.array(latvecs, dtype=float)
            self._dftbpluslib.dftbp_set_coords_and_lattice_vecs(
                self._dftb_handler, coords, latvecs)
        else:
            self._dftbpluslib.dftbp_set_coords(self._dftb_handler, coords)


    def set_external_potential(self, extpot, extpotgrad=None):
        '''Sets up an external potential.

        Args:
            extpot (2darray):     External potential at the position of each
                                  atom. Shape: [natom]. (in atomic units)
            extpotgrad (2darray): Gradient of the external potential at each
                                  atom. Shape: [natom, 3]. (in atomic units)
                                  This parameter is optional, you can pass None
                                  if you did not ask DFTB+ to calculate forces.
        '''

        if extpotgrad is not None:
            extpotgrad = np.ctypeslib.as_ctypes(extpotgrad)

        self._dftbpluslib.dftbp_set_external_potential(
            self._dftb_handler, extpot, extpotgrad)


    def get_nr_atoms(self):
        '''Queries the number of atoms.

        Returns:
            self._natoms (int): number of atoms
        '''

        return self._natoms


    def get_energy(self):
        '''Performs the energy (Mermin free) calculation
           and queries the energy of the current geometry.

        Returns:
            energy[0] (float): calculated Mermin free energy
                               (in atomic units)
        '''

        energy = np.empty(1)

        self._dftbpluslib.dftbp_get_energy(self._dftb_handler, energy)

        return energy[0]


    def get_gradients(self):
        '''Performs the calculation of the atomic gradients and
           queries the gradients of the current geometry.

        Returns:
            gradients (2darray): calculated gradients (in atomic units)
        '''

        gradients = np.empty((self._natoms, 3))

        self._dftbpluslib.dftbp_get_gradients(self._dftb_handler,
                                              gradients)

        return gradients


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

        grosschg = np.empty(self._natoms)

        self._dftbpluslib.dftbp_get_gross_charges(self._dftb_handler,
                                                  grosschg)

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
