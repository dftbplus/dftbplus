#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#
#
'''Convolves band.out-like files into DOS/PDOS via broadening functions'''

from __future__ import print_function
import argparse
import math
import numpy as np
import numpy.polynomial.hermite as  H
from dptools.bandout import BandOut
from dptools.scripts.common import ScriptError

USAGE = '''
Reads the band structure information stored in a file INPUT created by
DFTB+ (usually band.out for DOS and region*.out for PDOS) and
convolves the eigenlevels with a broadening function to produce nice
DOS/PDOS curves. The result is stored in the file OUTPUT in
NXY-format.

IMPORTANT: If you create PDOS, you have to weight the levels by the
'occupation' numbers (using option '-w'), otherwise you will obtain the total
DOS.

For spin unpolarized calculations, the output contains one Y-column only. For
spin polarized calculations the first column contains the total value, while
the further Y-columns contain the values for each spin-channel separately.
'''

DEFAULT_GRID_SEPARATION = 0.01

GAUSS_BROADENING = "gauss"
FERMI_BROADENING = "fermi"
MP_BROADENING = "mp"
BROADENING_FUNCTIONS = [GAUSS_BROADENING, FERMI_BROADENING, MP_BROADENING]

DEFAULT_RANGES = {
    GAUSS_BROADENING: 3.5,
    FERMI_BROADENING: 7.0,
    MP_BROADENING: 3.5,
}

DEFAULT_WIDTHS = {
    GAUSS_BROADENING: 0.1,
    FERMI_BROADENING: 0.1,
    MP_BROADENING: 0.1,
}


def main(cmdlineargs=None):
    '''Main driver for dp_dos.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    args, broadening, sigma, sigmarange, infile, outfile = \
    parse_arguments(cmdlineargs)
    dp_dos(args, broadening, sigma, sigmarange, infile, outfile)


def parse_arguments(cmdlineargs=None):
    '''Parses command line arguments.

    Args:
        cmdlineargs: List of command line arguments. When None, arguments in
            sys.argv are parsed (Default: None).
    '''
    parser = argparse.ArgumentParser(description=USAGE)

    msg = "grid separation (default: {0:.2f})".format(DEFAULT_GRID_SEPARATION)
    parser.add_argument("-g", "--gridres", type=float, dest="gridres",
                        default=DEFAULT_GRID_SEPARATION, help=msg)

    msg = "create pdos or occupation weighted dos"
    parser.add_argument("-w", "--weight-occ", action="store_true",
                        dest="occweight", default=False, help=msg)

    msg = "broadening function (default: gauss)"
    parser.add_argument("-f", "--broadening-function", dest="broadening",
                        choices=BROADENING_FUNCTIONS, default=GAUSS_BROADENING,
                        help=msg)

    msg = "order of mp function"
    parser.add_argument("-o", "--mporder", type=int, dest="mporder", help=msg)

    msg = "broadening width sigma (default: gauss {:.2f}, fermi {:.2f}, "\
          "mp {:.2f})".format(DEFAULT_WIDTHS[GAUSS_BROADENING],
                             DEFAULT_WIDTHS[FERMI_BROADENING],
                             DEFAULT_WIDTHS[MP_BROADENING])
    parser.add_argument("-b", "--broadening-width", type=float, metavar="WIDTH",
                        dest="broadwidth", help=msg, default=-1.0)

    msg = "number of sigmas after which the broadening function is considered "\
          "to be zero (default: gauss {0:.2f}, fermi {1:.2f}, mp {1:.2f})"\
           .format(DEFAULT_RANGES[GAUSS_BROADENING],
                   DEFAULT_RANGES[FERMI_BROADENING],
                   DEFAULT_RANGES[MP_BROADENING])
    parser.add_argument("-s", "--sigma-range", type=float, metavar="RANGE",
                        dest="broadrange", default=-1.0, help=msg)

    msg = "turn on verbose operation"
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                        default=False, help=msg)

    msg = "input file name"
    parser.add_argument("infile", metavar="INPUT", help=msg)

    msg = "output file name"
    parser.add_argument("outfile", metavar="OUTPUT", help=msg)

    args = parser.parse_args(cmdlineargs)

    if args.mporder:
        if args.broadening != 'mp':
            msg = '--mporder can only be set when --broadening-function=mp.'
            raise ScriptError(msg)
        if args.mporder < 1:
            raise ScriptError(msg)

    infile = args.infile
    outfile = args.outfile
    broadening = args.broadening

    if args.broadwidth < 0.0:
        sigma = DEFAULT_WIDTHS[broadening]
    else:
        sigma = args.broadwidth

    if args.broadrange < 0.0:
        sigmarange = DEFAULT_RANGES[broadening]
    else:
        sigmarange = args.broadrange

    return args, broadening, sigma, sigmarange, infile, outfile


def dp_dos(args, broadening, sigma, sigmarange, infile, outfile):
    '''convolves the eigenlevels with a broadening function to produce nice
       DOS/PDOS curves.

    Args:
        args: Containing the obtained parsed arguments.
        broadening: Specified broadening width sigma.
        sigma: Broadening width
        sigmarange: number of sigmas after which the broadening function is
                    considered to be zero
        infile: File containing the DFTB+ band structure information.
        outfile: Output file name
    '''
    try:
        bandout = BandOut.fromfile(infile)
    except OSError:
        raise ScriptError('You must enter a valid path to the input file.')
    eigvals = bandout.eigvalarray
    if not args.occweight:
        eigvals[:, :, :, 1] = 1.0

    if broadening == FERMI_BROADENING:
        aa = 0.5 / sigma
        bb = 1.0 / sigma
        broadening_function = lambda x: 1.0 / (1.0 + np.cosh(bb * x))
    elif broadening == GAUSS_BROADENING:
        aa = 1.0 / (sigma * np.sqrt(np.pi))
        bb = -1.0 / sigma**2
        broadening_function = lambda x: np.exp(bb * x * x)
    else: # Methfessel-Paxton
        aa = 1.0 / sigma
        bb = 1.0 / sigma
        order = args.mporder + 1 # order = 0 => Gaussian
        coefs = np.zeros(2*order, dtype=float)
        for ii in range(order):
            coefs[2*ii] = ((-1)**ii) / (math.factorial(ii) * (4**ii) \
                                        * math.sqrt(math.pi))
        broadening_function = lambda x: H.hermval(bb * x, coefs) \
                              * np.exp(-x * x * bb * bb)

    dsigma = sigmarange * sigma
    gridres = args.gridres

    minval = np.min(eigvals[:, :, :, 0]) - dsigma
    maxval = np.max(eigvals[:, :, :, 0]) + dsigma

    if args.verbose:
        if broadening == GAUSS_BROADENING:
            print("Gaussian broadening")
        elif broadening == FERMI_BROADENING:
            print("Fermi-function compatible broadening")
        else:
            print("Methfessel-Paxton compatible broadening")
        print("Broadening by {0:.2f} eV".format(sigma))
        if args.occweight:
            print("Weighting DOS by second column of data")
        print("Plotting from {0:.2f} eV to {1:.2f} eV".format(minval, maxval))
        print("Sigma window is : {0:.2f} eV with a grid resolution of {1:.2f}"
              .format(dsigma, gridres))

    # First and last grid points, x-grid
    xmin = np.floor((minval) / gridres) * gridres
    xmax = np.ceil((maxval) / gridres) * gridres
    xvals = np.arange(xmin, xmax + gridres, gridres, dtype=float)
    nval = len(xvals)

    # Empty container for y-values
    yvals = np.zeros((bandout.nspin, nval), dtype=float)

    # Calculate broadened values around each state on the grid.
    for ispin in range(bandout.nspin):
        for ik in range(bandout.nkpt):
            prefac = aa * bandout.kweights[ispin, ik]
            for eigval, occ in eigvals[ispin, ik]:
                # Grid points for the current curve (first, last)
                ilower = int(np.floor((eigval - dsigma - xmin) / gridres))
                iupper = int(np.ceil((eigval + dsigma - xmin) / gridres))
                dx = eigval - xvals[ilower:iupper+1]
                yvals[ispin, ilower:iupper+1] += (
                    (prefac * occ) * broadening_function(dx))
    # Write resulting DOS
    fp = open(outfile, "w")
    if bandout.nspin == 1:
        for xx, yy in zip(xvals, yvals[0]):
            fp.write("{0:18.10E} {1:18.10E}\n".format(xx, yy))
    else:
        ytotal = np.sum(yvals, axis=0)
        formstr0 = "{0:18.10E} "
        tmp = ["{" + "{0:d}".format(ii) + ":18.10E}" for ii in
               range(1, bandout.nspin + 2)]
        formstr = formstr0 + " ".join(tmp) + "\n"
        nvals = len(xvals)
        for ii in range(nvals):
            fp.write(formstr.format(xvals[ii], ytotal[ii], *yvals[:, ii]))
    fp.close()
