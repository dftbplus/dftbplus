#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
"""
Helper routines for finite difference scripts
"""

import numpy as np
import numpy.linalg as la
import shutil
import tempfile


BOHR__AA = 0.529177249
AA__BOHR = 1.0 / BOHR__AA


def readgen(fname):
    """Reads in the content of a gen file."""

    fp = open(fname, "r")
    line = fp.readline()
    words = line.split()
    natom = int(words[0])
    periodic = (words[1] == 'S' or words[1] == 's')
    fractional = (words[1] == 'F' or words[1] == 'f')
    periodic = (periodic or fractional)
    line = fp.readline()
    specienames = line.split()
    coords = np.empty((natom, 3), dtype=float)
    species = np.empty(natom, dtype=int)
    for ii in range(natom):
        line = fp.readline()
        words = line.split()
        species[ii] = int(words[1]) - 1
        coords[ii] = (float(words[2]), float(words[3]), float(words[4]))
    if periodic:
        line = fp.readline()
        origin = np.array([float(s) for s in line.split()], dtype=float)
        latvecs = np.empty((3, 3), dtype=float)
        for ii in range(3):
            line = fp.readline()
            latvecs[ii] = [float(s) for s in line.split()]
        if fractional:
            coords = frac2cart(latvecs, coords)
        origin *= AA__BOHR
        latvecs *= AA__BOHR
    else:
        origin = None
        latvecs = None
    coords *= AA__BOHR
    return specienames, species, coords, origin, latvecs


def writegen(fname, data):
    """Writes the geometry as gen file."""

    fp = open(fname, "w")
    specienames, species, coords, origin, latvecs = data
    fp.write("%5d %s\n" % (len(coords), latvecs is None and "C" or "S"))
    fp.write(("%2s "*len(specienames) + "\n") % tuple(specienames))
    coords = coords * BOHR__AA
    for ii in range(len(coords)):
        fp.write("%5d %5d %23.15E %23.15E %23.15E\n"
                 % (ii + 1, species[ii] + 1, coords[ii, 0], coords[ii, 1],
                    coords[ii, 2]))
    if latvecs is not None:
        origin = origin * BOHR__AA
        latvecs = latvecs * BOHR__AA
        fp.write("%23.15E %23.15E %23.15E\n" % tuple(origin))
        for ii in range(3):
            fp.write("%23.15E %23.15E %23.15E\n" % tuple(latvecs[ii]))
    fp.close()


def cart2frac(latvecs, coords):
    "Converts cartesian coordinates to fractional coordinates."

    invlatvecs = np.empty((3, 3), dtype=float)
    invlatvecs = np.transpose(latvecs)
    newcoords = np.array(coords)
    invlatvecs = la.inv(invlatvecs)
    for iat, atcoords in enumerate(coords):
        newcoords[iat] = np.dot(invlatvecs, atcoords)
    return newcoords


def frac2cart(latvecs, coords):
    """Converts fractional coordinates to cartesian ones."""

    newcoords = np.array(coords)
    for iat, atcoords in enumerate(coords):
        newcoords[iat] = np.dot(np.transpose(latvecs), atcoords)
    return newcoords


def stress2latderivs(stress, latvecs):
    """Converts stress to lattice derivatives."""

    invlatvecs = la.inv(latvecs)
    volume = la.det(latvecs)
    latderivs = -volume * np.transpose(np.dot(stress, invlatvecs))
    return latderivs


def create_temporary_copy(src_file_name):
    """
    Copies the source file into a temporary file.
    Returns a _TemporaryFileWrapper, whose destructor deletes the temp file
    (i.e. the temp file is deleted when the object goes out of scope).
    """
    tf = tempfile.NamedTemporaryFile()
    shutil.copy2(src_file_name, tf.name)
    return tf


def exists(filename):
    """
    Check for existence of named file
    """
    try:
        f = open(filename)
        f.close()
        return True
    except OSError:
        return False
