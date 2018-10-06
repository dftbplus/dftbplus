#!/usr/bin/env python
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
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

USAGE = '''Converts the given INPUT file in DFTB+ GEN format to CIF. Per default,
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
    infile, options = parse_cmdline_args(cmdlineargs)
    gen2cif(infile, options)

def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = argparse.ArgumentParser(description=USAGE, usage='%(prog)s [options] INPUT')
    parser.add_argument("-o", "--output", action="store", dest="output",
                        help="override the name of the output file (use '-' for "
                        "standard out")
    parser.add_argument("-c", "--cellsize", action="store", dest="cellsize",
                        type=float, default=100.0, help="lattice constant "
                        "for the simple cubic cell created for non-periodic "
                        " geometries (default: 100)")
    options, args = parser.parse_known_args(cmdlineargs)

    if len(args) != 1:
        raise ScriptError('You must specify exactly one argument (input file).')
    infile = args[0]

    return infile, options

def gen2cif(infile, options):
    '''Converts the given INPUT file in DFTB+ GEN format to CIF format.

    Args:
        infile: File containing the gen-formatted geometry.
        options: Options (e.g. as returned by the command line parser).
    '''
    try:
        gen = Gen.fromfile(infile)
    except OSError:
        raise ScriptError('You must enter a valid path to the input file.')

    geometry = gen.geometry
    if not geometry.periodic:
        geometry.setlattice(options.cellsize * np.eye(3))
    cif = Cif(gen.geometry)

    if options.output:
        if options.output == "-":
            outfile = sys.stdout
        else:
            outfile = options.output
    else:
        if infile.endswith(".gen"):
            outfile = infile[:-4] + ".cif"
        else:
            outfile = infile + ".cif"
    cif.tofile(outfile)
