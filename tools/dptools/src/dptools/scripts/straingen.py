#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2017  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Applies strain to a DFTB+ gen file.'''

import sys
import optparse
import numpy as np
from dptools.gen import Gen
from dptools.scripts.common import ScriptError

USAGE = """usage: %prog [options] INPUT

Strains the geometry found in INPUT, writing the resulting geometries
to standard output."""

# Voight convention for 1 index to 2 index for strain tensors
VOIGHT = [[0, 0], [1, 1], [2, 2], [1, 2], [0, 2], [0, 1]]
# Labels for the types of strain
LABELS = {'xx': (0, ), 'yy': (1, ), 'zz': (2, ), 'yz': (3, ), 'xz': (4, ),
          'xy': (5, ), 'i': (0, 1, 2)}


def main(cmdlineargs=None):
    '''Main driver for straingen.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed. (Default: None)
    '''
    infile, options = parse_cmdline_args(cmdlineargs)
    straingen(infile, options)

def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed. (Default: None)
    '''
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-o", "--output", action="store", dest="output", default='-',
                      help="override the name of the output file (use '-' for "
                      "standard out")
    parser.add_option("-s", "--strain", action="store", dest="strain",
                      type=float, default=0.0, help="positive or negative "
                      "percentage strain for the geometries (default: 0)")
    parser.add_option("-c", "--component", action="store", dest="component",
                      type=str, default='I', help="strain type to apply "
                      "posible values being xx, yy, zz, xz, xz, yz or I for "
                      "isotropic (default value: I)")

    options, args = parser.parse_args(cmdlineargs)

    if options.component.lower() not in LABELS:
        msg = "Invalid strain component '" + options.component + "'"
        raise ScriptError(msg)

    if len(args) != 1:
        raise ScriptError("You must specify exactly one argument (input file).")
    infile = args[0]

    return infile, options

def straingen(infile, options):
    '''Strains a geometry from a gen file.

    Args:
        infile: File containing the gen-formatted geometry
        options: Options (e.g. as returned by the command line parser)
    '''

    gen = Gen.fromfile(infile)
    geometry = gen.geometry

    strain = np.zeros((3, 3), dtype=float)
    for jj in range(3):
        strain[jj][jj] = 1.0

    components = LABELS[options.component.lower()]

    for ii in components:
        strain[VOIGHT[ii][0]][VOIGHT[ii][1]] += 0.005*options.strain
        strain[VOIGHT[ii][1]][VOIGHT[ii][0]] += 0.005*options.strain

    if geometry.latvecs is not None:
        geometry.latvecs = np.dot(geometry.latvecs, strain)

    geometry.coords = np.dot(geometry.coords, strain)

    if options.output:
        if options.output == "-":
            outfile = sys.stdout
        else:
            outfile = options.output
    else:
        if infile.endswith(".gen"):
            outfile = infile
        else:
            outfile = infile + ".gen"

    gen = Gen(geometry, fractional=gen.fractional)
    gen.tofile(outfile)
