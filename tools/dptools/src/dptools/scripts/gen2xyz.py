#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Convert DFTB+ gen format to XYZ.'''

import sys
import argparse
from dptools.gen import Gen
from dptools.xyz import Xyz
from dptools.scripts.common import ScriptError

USAGE = '''
Converts the given INPUT file in DFTB+ GEN format to XYZ. Per default,
if the filename INPUT is of the form PREFIX.gen the result is stored in PREFIX.xyz,
otherwise in INPUT.xyz. You can additionally store lattice vectors of the GEN
file in a separate file.
'''

def main(cmdlineargs=None):
    '''Main driver routine for gen2xyz.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    args = parse_cmdline_args(cmdlineargs)
    gen2xyz(args)

def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = argparse.ArgumentParser(description=USAGE)
    msg = "store lattice vectors in an external file"
    parser.add_argument("-l", "--lattice-file", action="store", dest="lattfile",
                        help=msg)
    msg = "override the name of the output file (use '-' for standard output)"
    parser.add_argument("-o", "--output", action="store", dest="output",
                        help=msg)
    msg = "comment for the second line of the xyz-file"
    parser.add_argument("-c", "--comment", action="store", dest="comment",
                        default="", help=msg)
    msg = "input file name"
    parser.add_argument("infile", metavar="INPUT", help=msg)

    args = parser.parse_args(cmdlineargs)

    return args

def gen2xyz(args):
    '''Converts the given INPUT file in DFTB+ GEN format to XYZ format.

    Args:
        args: Namespace of command line arguments
    '''
    infile = args.infile
    try:
        gen = Gen.fromfile(infile)
    except OSError:
        raise ScriptError('You must enter a valid path to the input file.')
    xyz = Xyz(gen.geometry, args.comment)
    if args.output:
        if args.output == "-":
            outfile = sys.stdout
        else:
            outfile = args.output
    else:
        if infile.endswith(".gen"):
            outfile = infile[:-4] + ".xyz"
        else:
            outfile = infile + ".xyz"
    xyz.tofile(outfile)

    if gen.geometry.periodic and args.lattfile:
        fp = open(args.lattfile, "w")
        for vec in gen.geometry.latvecs:
            fp.write("{0:18.10E} {1:18.10E} {2:18.10E}\n".format(*vec))
        fp.close()
