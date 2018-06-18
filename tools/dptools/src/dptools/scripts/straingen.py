#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
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

USAGE = """usage: %prog [options] INPUT STRAIN

Strains the geometry found in INPUT, writing the resulting geometries
to standard output. Possible values for the strain type to apply are
xx, yy, zz, xz, xz, yz or I for isotropic. STRAIN is specified as
positive or negative percentage for the geometries.

Note: in case of negative strain you have to separately prefix the
argument block by '--', e.g. '%prog -c I -- geo.gen -5'"""

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
    infile, strain, options = parse_cmdline_args(cmdlineargs)
    straingen(infile, strain, options)


def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed. (Default: None)
    '''
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option("-o", "--output", action="store", dest="output",
                      default='-', help="override the name of the output file "
                      "(use '-' for standard out")
    parser.add_option("-c", "--component", action="store", dest="component",
                      type=str, default='I', help="strain type to apply "
                      "posible values being xx, yy, zz, xz, xz, yz or I for "
                      "isotropic (default value: I)")
    options, args = parser.parse_args(cmdlineargs)

    if options.component.lower() not in LABELS:
        msg = "Invalid strain component '" + options.component + "'"
        raise ScriptError(msg)

    if len(args) != 2:
        raise ScriptError("You must specify exactly two arguments "
                          "(input file, strain).")
    infile = args[0]
    try:
        strain = float(args[1])
    except ValueError:
        raise ScriptError("Invalid strain value '" + args[1] + "'")

    if options.component.lower() in ('xx', 'yy', 'zz', 'i') and strain <= -100:
        raise ScriptError("Compressive strain cannot exceed 100%")

    return infile, strain, options


def straingen(infile, strain, options):
    '''Strains a geometry from a gen file.

    Args:
        infile: File containing the gen-formatted geometry
        strain: Strain to apply
        options: Options (e.g. as returned by the command line parser)
    '''

    gen = Gen.fromfile(infile)
    geometry = gen.geometry

    strainmtx = np.zeros((3, 3), dtype=float)
    for jj in range(3):
        strainmtx[jj][jj] = 1.0

    components = LABELS[options.component.lower()]

    for ii in components:
        strainmtx[VOIGHT[ii][0]][VOIGHT[ii][1]] += 0.005 * strain
        strainmtx[VOIGHT[ii][1]][VOIGHT[ii][0]] += 0.005 * strain

    if geometry.latvecs is not None:
        geometry.latvecs = np.dot(geometry.latvecs, strainmtx)

    geometry.coords = np.dot(geometry.coords, strainmtx)

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
