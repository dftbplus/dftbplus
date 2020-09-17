#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Converts DFTB+ gen format to CIF.'''

import sys
import argparse
import numpy as np
from dptools.gen import Gen
from dptools.cif import Cif
from dptools.scripts.common import ScriptError

USAGE = '''
Converts the given INPUT file in DFTB+ GEN format to CIF. Per default,
if the filename INPUT is of the form PREFIX.gen the result is stored in PREFIX.cif,
otherwise in INPUT.cif. Since the GEN format does not contain any symmetry
information, the symmetry is set to P1 in the CIF file. If the GEN format
contains a non-periodic geometry, the lattice in the CIF format is set to
simple cubic.
'''

def main(cmdlineargs=None):
    '''Main driver routine of gen2cif.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    args = parse_cmdline_args(cmdlineargs)
    gen2cif(args)

def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = argparse.ArgumentParser(description=USAGE)
    msg = "override the name of the output file (use '-' for standard out)"
    parser.add_argument("-o", "--output", action="store", dest="output",
                        help=msg)
    msg = "lattice constant for the simple cubic cell created for non-periodic"\
          " geometries (default: 100)"
    parser.add_argument("-c", "--cellsize", action="store", dest="cellsize",
                        type=float, default=100.0, help=msg)
    msg = "input file name"
    parser.add_argument("infile", metavar="INPUT", help=msg)

    args = parser.parse_args(cmdlineargs)

    return args

def gen2cif(args):
    '''Converts the given INPUT file in DFTB+ GEN format to CIF format.

    Args:
        args: Namespace of command line arguments
    '''
    infile = args.infile
    try:
        gen = Gen.fromfile(infile)
    except OSError:
        raise ScriptError('You must enter a valid path to the input file.')

    geometry = gen.geometry
    if not geometry.periodic:
        geometry.setlattice(args.cellsize * np.eye(3))
    cif = Cif(gen.geometry)

    if args.output:
        if args.output == "-":
            outfile = sys.stdout
        else:
            outfile = args.output
    else:
        if infile.endswith(".gen"):
            outfile = infile[:-4] + ".cif"
        else:
            outfile = infile + ".cif"
    cif.tofile(outfile)
