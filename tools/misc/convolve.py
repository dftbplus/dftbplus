#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

"""
Broadens the eigenvalues with Gaussians of given width
"""

from numpy import size, arange
from math import pi, sqrt, exp
import re
import argparse
import sys


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="infile", required=True,
                    help="input file name")
parser.add_argument("-o", "--output", dest="outfile", required=True,
                    help="output file name")
parser.add_argument("-b", "--broaden", dest="sigma", type=float, default=0.1,
                    help="broadening width")
parser.add_argument("-u", "--unnorm", dest="norm", action="store_false",
                    help="don't normalize by number of spins/kpoints")
parser.add_argument("-w", "--weight", dest="weight", action="store_true",
                    help="weight data present")

args = parser.parse_args()

outfile = args.outfile
infile = args.infile
sigma = args.sigma
norm = args.norm
weight = args.weight

# returns numeric matches as a number
def numGrep(pattern,fileObj):
    r=[]
    for line in fileObj:
        m = re.match(pattern,line)
        if m:
            r.append(m.group(1,2))
    return r

gridResltn = 3
figures = "%i" % gridResltn
sGrid = "%."+figures+"f"
gridResltn = 1.0/(10**gridResltn)

print("Gaussian broadening of levels, width %.2f" %(sigma))

a = 1.0 / (sigma * sqrt(2.0*pi))
b = - 0.5 / (sigma * sigma)

states = {}

fp = open(infile, "r")
eigenstates = numGrep(r"^ *([+-]?\d+\.\d+) +([+-]?\d+\.\d+e?-?\d*)",fp)
fp.close()

if norm:
    fp = open(infile, "r")
    nKPTS = numGrep(r"^ *KPT *(\d+) *(.*)",fp)
    scale = 2.0 / size(nKPTS)
    fp.close()
else:
    scale = 1.0

for n in eigenstates:
    x = float(n[0])
    C = float(n[1])
    lower = float(sGrid % (x-3.5*sigma))
    upper = float(sGrid % (x+3.5*sigma))
    if weight:
        C *= a
    else:
        C = a
    C *= scale
    for m in arange(lower,upper,gridResltn):
        val = float(sGrid % (m))
        if val in states:
            states[val] += C * exp(b * (x-m)*(x-m))
        else:
            states[val] = C * exp(b * (x-m)*(x-m))

fp = open(outfile, "w")
lower = min(states.keys())
upper = max(states.keys())
for m in arange(lower,upper,gridResltn):
    val = float(sGrid % (m))
    if val in states:
        fp.write("%s %s\n" % (m,states[val]))
    else:
        fp.write("%s 0.0\n" % (m))
fp.close()
