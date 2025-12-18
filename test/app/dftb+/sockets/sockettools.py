#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

import numpy as np

a0 = 0.529177249

def frac2cart(latvecs,coords):
    newcoords = np.array(coords)
    for iAt in range(len(coords)):
        newcoords[iAt] = np.dot(np.transpose(latvecs),coords[iAt])
    return (newcoords)


def readgen(fname):
    fp = open(fname, "r")
    line = fp.readline()
    words = line.split()
    nAtom = int(words[0])
    periodic = (words[1] == 'S' or words[1] == 's')
    fractional = (words[1] == 'F' or words[1] == 'f')
    periodic = (periodic or fractional)
    line = fp.readline()
    specienames = line.split()
    coords = np.empty((nAtom, 3), dtype=float)
    species = np.empty(nAtom, dtype=int)
    for ii in range(nAtom):
        line = fp.readline()
        words = line.split()
        species[ii] = int(words[1]) - 1
        coords[ii] = (float(words[2]), float(words[3]), float(words[4]))
    if periodic:
        line = fp.readline()
        origin = np.array([ float(s) for s in line.split() ], dtype=float)
        latvecs = np.empty((3, 3), dtype=float)
        for ii in range(3):
            line = fp.readline()
            latvecs[ii] = [ float(s) for s in line.split() ]
        if fractional :
            coords = frac2cart(latvecs,coords)
    else:
        origin = None
        latvecs = None
    return specienames, species, coords, origin, latvecs


def receive_all(connection, msglen):
    received = b''
    nn = 0
    while nn < msglen:
        fragment = connection.recv(msglen - nn)
        if fragment == '':
            raise RuntimeError("socket connection broken")
        nn += len(fragment)
        received += fragment
    return received
