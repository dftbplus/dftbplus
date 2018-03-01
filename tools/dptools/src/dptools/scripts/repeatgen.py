#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Repeats a geometry along supercell vectors'''

import sys
import optparse
import numpy as np
from dptools.gen import Gen
from dptools.geometry import Geometry
from dptools.scripts.common import ScriptError

USAGE = '''usage: %prog [options] INPUT N1 N2 N3

Repeats the geometry found in INPUT along each supercell vector N1, N2
and N3 times, respectively and writes the resulting geometry to
standard output'''


def main(cmdlineargs=None):
    '''Repeatgen main driver.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    infile, repeats, options = parse_cmdline_args(cmdlineargs)
    repeatgen(infile, repeats, options)


def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = optparse.OptionParser(usage=USAGE)
    msg = 'file containing lattice vectors (overrides lattice vectors'\
          ' in the geometry file)'
    parser.add_option(
        '-l', '--lattice-file', action='store', help=msg, dest='latticefile')
    msg = 'output file to store the resulting geometry'
    parser.add_option('-o', '--output', action='store', default='-', help=msg)
    options, args = parser.parse_args(cmdlineargs)

    if len(args) != 4:
        raise ScriptError('Incorrect number of arguments')
    infile = args[0]
    reps = []
    for repstr in args[1:4]:
        try:
            reps.append(int(repstr))
        except ValueError:
            msg = "Invalid repetition number '" + repstr + "'."
            raise ScriptError(msg)
    if not (reps[0] > 0 and reps[1] > 0 and reps[2] > 0):
        raise ScriptError('Repetition numbers must be greater than zero')

    return infile, (reps[0], reps[1], reps[2]), options


def repeatgen(infile, repeats, options):
    '''Repeats geometry from gen files.

    Args:
        infile: File containing the gen-formatted geometry.
        repeats: (n1, n2, n3) integer tuple containing the repetitions along
            each lattice vector.
        options: Options (e.g. as returned by the command line parser).
    '''

    gen = Gen.fromfile(infile)
    geo = gen.geometry

    latvecs = geo.latvecs
    if options.latticefile:
        latvecs = np.fromfile(options.latticefile, sep=' ')
        if len(latvecs) != 9:
            msg = ('Invalid number of lattice vector components in '
                   + options.latticefile)
            raise ScriptError(msg)
        latvecs.shape = (3, 3)

    if latvecs is None:
        msg = 'No lattice vectors found (neither in gen nor in external file)'
        raise ScriptError(msg)

    newgeo = _repeatgeo(geo, latvecs, repeats)
    newgen = Gen(newgeo, gen.fractional)
    outfile = options.output
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
