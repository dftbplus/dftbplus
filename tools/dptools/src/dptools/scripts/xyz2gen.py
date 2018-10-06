#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Converts XYZ to DFTB+ gen format'''

import sys
import argparse
import numpy as np
from dptools.gen import Gen
from dptools.xyz import Xyz
from dptools.scripts.common import ScriptError

USAGE = '''Converts the given INPUT file in XYZ format to DFTB+ GEN format. Per default, if
the filename INPUT is of the form PREFIX.xyz the result is stored in PREFIX.gen,
otherwise in INPUT.gen. You can additionally specify a file with lattice
vectors to create a periodic structure in the GEN file.'''


def main(cmdlineargs=None):
    '''Main driver routine of xyz2gen.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    infile, options = parse_cmdline_args(cmdlineargs)
    xyz2gen(infile, options)

def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = argparse.ArgumentParser(description=USAGE, usage='%(prog)s [options] INPUT')
    parser.add_argument("-l", "--lattice-file", action="store", dest="lattfile",
                        help="read lattice vectors from an external file")
    parser.add_argument("-o", "--output", action="store", dest="output",
                        help="override the name of the output file (use '-' for "
                        "standard out")
    parser.add_argument("-f", "--fractional", action="store_true",
                        dest="fractional", default=False,
                        help="store coordinate in fractional format instead of "
                        "absolute coordinates")
    options, args = parser.parse_known_args(cmdlineargs)

    if len(args) != 1:
        raise ScriptError('You must specify exactly one argument (input file).')
    infile = args[0]

    return infile, options

def xyz2gen(infile, options):
    '''Converts the given INPUT file in XYZ format to DFTB+ GEN format.

    Args:
        infile: File containing the xyz-formatted geometry.
        options: Options (e.g. as returned by the command line parser).
    '''
    try:
        xyz = Xyz.fromfile(infile)
    except OSError:
        raise ScriptError('You must enter a valid path to the input file.')
    geo = xyz.geometry
    if options.lattfile:
        fp = open(options.lattfile, "r")
        tmp = np.fromfile(fp, count=9, dtype=float, sep=" ")
        latvecs = tmp.reshape((3, 3))
        fp.close()
        geo.setlattice(latvecs)
    gen = Gen(geo, fractional=options.fractional)

    if options.output:
        if options.output == "-":
            outfile = sys.stdout
        else:
            outfile = options.output
    else:
        if infile.endswith(".xyz"):
            outfile = infile[:-4] + ".gen"
        else:
            outfile = infile + ".gen"
    gen.tofile(outfile)
