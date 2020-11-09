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

USAGE = '''
Converts the given INPUT file in XYZ format to DFTB+ GEN format. Per default, if
the filename INPUT is of the form PREFIX.xyz the result is stored in PREFIX.gen,
otherwise in INPUT.gen. You can additionally specify a file with lattice
vectors to create a periodic structure in the GEN file.
'''


def main(cmdlineargs=None):
    '''Main driver routine of xyz2gen.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    args = parse_cmdline_args(cmdlineargs)
    xyz2gen(args)

def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = argparse.ArgumentParser(description=USAGE)
    parser.add_argument("-l", "--lattice-file", action="store", dest="lattfile",
                        help="read lattice vectors from an external file")
    parser.add_argument("-o", "--output", action="store", dest="output",
                        help="override the name of the output file (use '-' for"
                        " standard out)")
    parser.add_argument("-f", "--fractional", action="store_true",
                        dest="fractional", default=False,
                        help="store coordinate in fractional format instead of "
                        "absolute coordinates")
    msg = "input file name"
    parser.add_argument("infile", metavar="INPUT", help=msg)

    args = parser.parse_args(cmdlineargs)

    return args

def xyz2gen(args):
    '''Converts the given INPUT file in XYZ format to DFTB+ GEN format.

    Args:
        args: Namespace of command line arguments
    '''
    infile = args.infile
    try:
        xyz = Xyz.fromfile(infile)
    except OSError:
        raise ScriptError('You must enter a valid path to the input file.')
    geo = xyz.geometry
    if args.lattfile:
        fp = open(args.lattfile, "r")
        tmp = np.fromfile(fp, count=9, dtype=float, sep=" ")
        latvecs = tmp.reshape((3, 3))
        fp.close()
        geo.setlattice(latvecs)
    gen = Gen(geo, fractional=args.fractional)

    if args.output:
        if args.output == "-":
            outfile = sys.stdout
        else:
            outfile = args.output
    else:
        if infile.endswith(".xyz"):
            outfile = infile[:-4] + ".gen"
        else:
            outfile = infile + ".gen"
    gen.tofile(outfile)
