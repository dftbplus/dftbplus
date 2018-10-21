#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#

'''Applies strain to a DFTB+ gen file.'''

import sys
import argparse
import numpy as np
from dptools.gen import Gen
from dptools.scripts.common import ScriptError

USAGE = '''
Strains the geometry found in INPUT, writing the resulting geometries
to standard output. Possible values for the strain type to apply are xx, yy, zz,
xz, xz, yz or I for isotropic. STRAIN is specified as positive or negative
percentage for the geometries.
'''

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
    args, strain = parse_cmdline_args(cmdlineargs)
    straingen(args, strain)


def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed. (Default: None)
    '''
    parser = argparse.ArgumentParser(description=USAGE)
    msg = "override the name of the output file (use '-' for standard out)"
    parser.add_argument("-o", "--output", action="store", dest="output",
                        default='-', help=msg)
    msg = "strain type to apply posible values being xx, yy, zz, xz, xz, yz or"\
          " I for isotropic (default value: I)"
    parser.add_argument("-c", "--component", action="store", dest="component",
                        type=str, default='I', help=msg)
    msg = "input file name"
    parser.add_argument("infile", metavar="INPUT", help=msg)
    msg = "strain to apply as positive or negative percentage"
    parser.add_argument("strain", type=float, metavar="STRAIN", help=msg)

    args = parser.parse_args(cmdlineargs)

    if args.component.lower() not in LABELS:
        msg = "Invalid strain component '" + args.component + "'"
        raise ScriptError(msg)

    strain = args.strain

    if args.component.lower() in ('xx', 'yy', 'zz', 'i') and strain <= -100:
        raise ScriptError("Compressive strain cannot exceed 100%")

    return args, strain


def straingen(args, strain):
    '''Strains a geometry from a gen file.

    Args:
        args: Namespace of command line arguments
        strain: Strain to apply
    '''
    infile = args.infile
    try:
        gen = Gen.fromfile(infile)
    except OSError:
        raise ScriptError('You must enter a valid path to the input file.')
    geometry = gen.geometry

    strainmtx = np.zeros((3, 3), dtype=float)
    for jj in range(3):
        strainmtx[jj][jj] = 1.0

    components = LABELS[args.component.lower()]

    for ii in components:
        strainmtx[VOIGHT[ii][0]][VOIGHT[ii][1]] += 0.005 * strain
        strainmtx[VOIGHT[ii][1]][VOIGHT[ii][0]] += 0.005 * strain

    if geometry.latvecs is not None:
        geometry.latvecs = np.dot(geometry.latvecs, strainmtx)

    geometry.coords = np.dot(geometry.coords, strainmtx)

    if args.output:
        if args.output == "-":
            outfile = sys.stdout
        else:
            outfile = args.output
    else:
        if infile.endswith(".gen"):
            outfile = infile
        else:
            outfile = infile + ".gen"

    gen = Gen(geometry, fractional=gen.fractional)
    gen.tofile(outfile)
