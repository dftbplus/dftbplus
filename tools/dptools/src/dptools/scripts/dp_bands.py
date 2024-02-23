#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2023  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Converts band.out to plottable data.'''

from argparse import ArgumentParser
import numpy as np
from dptools.bandout import BandOut
from dptools.scripts.common import ScriptError, find_auto_alignment

USAGE = '''
Reads the band structure information stored in the file INPUT created
by DFTB+ (usually band.out) and converts it to NXY-data (to be visualized by
gnuplot, xmgrace, etc.). The outputfile will be called OUTPREFIX_tot.dat and
contains the band structure for all spin channels. If spin channel separation is
specified, the output files OUTPREFIX_s1.dat, OUTPREFIX_s2.dat, etc. will be
created for each spin channel.
'''

def main(cmdlineargs=None):
    '''Main driver for dp_bands.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    args = parse_cmdline_args(cmdlineargs)
    dp_bands(args)


def parse_cmdline_args(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = ArgumentParser(description=USAGE)

    msg = "do not use the first column of the output to enumerate the k-points"
    parser.add_argument("-N", "--no-enumeration", action="store_false",
                        default=True, dest="enum", help=msg)

    msg = "create separate band structure for each spin channel"
    parser.add_argument("-s", "--separate-spins", action="store_true",
                        default=False, dest="spinsep", help=msg)

    msg = "Aligns an estimate of the level of half-occupation to the zero of" +\
        " energy. Integer occupations: Shifts the VBM to zero energy." +\
        " Fractional occupations: Shifts approximate Fermi level, i.e. the highest" +\
        " eigenvalue that is at least half occupied, to zero energy." +\
        " Note that this energy does not necessarily correspond to the" +\
        " real Fermi level. Please consult the appropriate output files" +\
        " (detailed.out|results.tag) to obtain the Fermi level of the" +\
        " calculation."
    parser.add_argument("-A", "--auto-align", action="store_true",
                        default=False, dest="autoalign", help=msg)

    msg = "manually defines zero energy [eV]"
    parser.add_argument("-a", "--align", type=float, default=None, help=msg)

    msg = "input file name"
    parser.add_argument("infile", metavar="INPUT", help=msg)

    msg = "output prefix"
    parser.add_argument("outprefix", metavar="OUTPREFIX", help=msg)

    args = parser.parse_args(cmdlineargs)

    # check for argument collisions
    if args.autoalign and args.align is not None:
        raise ScriptError(
            'Collision of incompatible options. Auto-alignment ' + \
            'and simultaneous manual shifting of eigenvalues is not supported.')

    return args


def dp_bands(args):
    '''Converts band structure information of DFTB+ output to NXY-format.

    Args:
        args: Namespace of command line arguments
    '''
    infile = args.infile
    outprefix = args.outprefix

    try:
        bandout = BandOut.fromfile(infile)
    except OSError:
        raise ScriptError('You must enter a valid path to the input file.') \
            from OSError

    if args.spinsep:
        # Create filename for each spin-channel
        fnames = ["{0}_s{1:d}.dat".format(outprefix, ispin)
                  for ispin in range(1, bandout.nspin + 1)]
        plotvals = bandout.eigvalarray[:, :, :, 0]
    else:
        # Put all eigenvalues in one channel, if no spin separation required
        eigvals = bandout.eigvalarray[:, :, :, 0]
        plotvals = np.empty((1, eigvals.shape[1],
                             eigvals.shape[2] * eigvals.shape[0]), dtype=float)
        for ispin in range(bandout.nspin):
            istart = ispin * bandout.nstate
            iend = (ispin + 1) * bandout.nstate
            plotvals[0, :, istart:iend] = eigvals[ispin, :, :]
        fnames = [outprefix + "_tot.dat",]

    if args.autoalign:
        plotvals -= find_auto_alignment(bandout)

    if args.align is not None:
        plotvals -= args.align

    if args.enum:
        formstr0 = "{0:d} "
        tmp = ["{" + str(ii) + ":18.10E}"
               for ii in range(1, plotvals.shape[2] + 1)]
        formstr = formstr0 + " ".join(tmp) + "\n"
    else:
        tmp = ["{" + str(ii) + ":18.10E}"
               for ii in range(plotvals.shape[2])]
        formstr = " ".join(tmp) + "\n"

    for iout, data in enumerate(plotvals):
        fp = open(fnames[iout], "w")
        for ik, kdata in enumerate(data):
            if args.enum:
                fp.write(formstr.format(ik + 1, *kdata))
            else:
                fp.write(formstr.format(*kdata))
        fp.close()
