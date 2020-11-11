#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Repeats a geometry along supercell vectors'''

import sys
import argparse
import numpy as np
from dptools.gen import Gen
from dptools.geometry import Geometry
from dptools.scripts.common import ScriptError

USAGE = '''
Repeats the geometry found in INPUT along each supercell vector N1,
N2 and N3 times, respectively and writes the resulting geometry to standard
output
'''


def main(cmdlineargs=None):
    '''Repeatgen main driver.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    args = parse_cmdline_args(cmdlineargs)
    repeatgen(args)


def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = argparse.ArgumentParser(description=USAGE)
    msg = 'file containing lattice vectors (overrides lattice vectors'\
          ' in the geometry file)'
    parser.add_argument(
        '-l', '--lattice-file', action='store', help=msg, dest='latticefile')
    msg = 'output file to store the resulting geometry'
    parser.add_argument('-o', '--output', action='store', default='-', help=msg)
    msg = 'create a repeat geometry for phonon bandstructure'
    parser.add_argument('-p', '--phonons', action='store_true', default=False,
                        help=msg)
    msg = 'input file name'
    parser.add_argument("infile", metavar="INPUT", help=msg)
    msg = 'repetition along the first lattice vector'
    parser.add_argument('n1', type=int, metavar="N1", help=msg)
    msg = 'repetition along the second lattice vector'
    parser.add_argument('n2', type=int, metavar="N2", help=msg)
    msg = 'repetition along the third lattice vector'
    parser.add_argument('n3', type=int, metavar="N3", help=msg)

    args = parser.parse_args(cmdlineargs)

    if not (args.n1 > 0 and args.n2 > 0 and args.n3 > 0):
        raise ScriptError('Repetition numbers must be greater than zero')
    if args.phonons:
        if (args.n1 % 2 ==0 or args.n2 % 2 == 0 or args.n3 % 2 == 0):
            raise ScriptError('Repetition numbers must be odd numbers')


    return args


def repeatgen(args):
    '''Repeats geometry from gen files.

    Args:
        args: Namespace of command line arguments
    '''
    infile = args.infile
    repeats = [args.n1, args.n2, args.n3]

    try:
        gen = Gen.fromfile(infile)
    except OSError:
        raise ScriptError('You must enter a valid path to the input file.')
    geo = gen.geometry

    latvecs = geo.latvecs
    if args.latticefile:
        latvecs = np.fromfile(args.latticefile, sep=' ')
        if len(latvecs) != 9:
            msg = ('Invalid number of lattice vector components in '
                   + args.latticefile)
            raise ScriptError(msg)
        latvecs.shape = (3, 3)

    if latvecs is None:
        msg = 'No lattice vectors found (neither in gen nor in external file)'
        raise ScriptError(msg)

    if args.phonons:
        newgeo = _repeatgeo2(geo, latvecs, repeats)
    else:
        newgeo = _repeatgeo(geo, latvecs, repeats)

    newgen = Gen(newgeo, gen.fractional)
    outfile = args.output
    if outfile == '-':
        outfile = sys.stdout
    newgen.tofile(outfile)


def _repeatgeo(geo, latvecs, repeats):
    '''Repeats geometry along given lattice vectors'''
    natoms = geo.natom
    coords = geo.coords
    inds = geo.indexes
    images = repeats[0] * repeats[1] * repeats[2]
    allcoords = np.empty((images * natoms, 3), dtype=float)
    allcoords[0:natoms, :] = coords
    ind = 0
    currepeats = np.zeros(3, dtype=int)
    for currepeats[0] in range(repeats[0]):
        for currepeats[1] in range(repeats[1]):
            for currepeats[2] in range(repeats[2]):
                shift = np.sum(latvecs * currepeats[:, np.newaxis], axis=0)
                shiftedcoords = coords + shift
                allcoords[ind * natoms : (ind + 1) * natoms, :] = shiftedcoords
                ind += 1
    allinds = np.tile(inds, images)
    repeats = np.array(repeats)
    if geo.periodic:
        newlatvecs = latvecs * repeats[:, np.newaxis]
    else:
        newlatvecs = None
    newgeo = Geometry(geo.specienames, allinds, allcoords, latvecs=newlatvecs,
                      origin=geo.origin)
    return newgeo

def _repeatgeo2(geo, latvecs, repeats):
    '''Repeats geometry along given lattice vectors for phonon calculations'''
    natoms = geo.natom
    coords = geo.coords
    inds = geo.indexes
    rep = np.array([(repeats[0]-1)/2, (repeats[1]-1)/2, (repeats[2]-1)/2])
    images = (repeats[0] * repeats[1] * repeats[2])
    allcoords = np.empty((images * natoms, 3), dtype=float)
    allcoords[0:natoms, :] = coords    
    ind = 1
    currepeats = np.zeros(3, dtype=int)
    for currepeats[2] in range(-rep[2],rep[2]+1):
        for currepeats[1] in range(-rep[1],rep[1]+1):
            for currepeats[0] in range(-rep[0],rep[0]+1):
                if (currepeats[0]==0 and currepeats[1]==0 and currepeats[2]==0):
                    continue
                shift = np.sum(latvecs * currepeats[:, np.newaxis], axis=0)
                shiftedcoords = coords + shift
                allcoords[ind * natoms : (ind + 1) * natoms, :] = shiftedcoords
                ind += 1
    allinds = np.tile(inds, images)
    repeats = np.array(repeats)
    if geo.periodic:
        newlatvecs = latvecs * repeats[:, np.newaxis]
    else:
        newlatvecs = None
    newgeo = Geometry(geo.specienames, allinds, allcoords, latvecs=newlatvecs,
                      origin=geo.origin)
    return newgeo

