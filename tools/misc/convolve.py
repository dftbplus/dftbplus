#!/usr/bin/env python3
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2020  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

from numpy import *
from math import *
import re
import getopt, sys

sigma = 0.1 # Generally looks OK
weight = False
outfile = None
infile = None
norm = True

def useage():
    print("--broaden -b broadening width")
    print("--help    -h this message")
    print("--input   -i input file name")
    print("--unnorm  -u don't normalize by number of spins/kpoints")
    print("--output  -o output file name")
    print("--weight  -w weight data present")
    sys.exit(2)

try:
    opts, args = getopt.getopt(sys.argv[1:], "o:i:b:hwu", ["output=", "input=", "broaden=", "help","weight","unnorm"])
except getopt.GetoptError as err:
    print(str(err)) # will print something like "option -a not recognized"
    useage()    
output = None
input = None
for options, argument in opts:
    if options in ("-h", "--help"):
        useage()
    elif options in ("-o", "--output"):
        outfile = argument
    elif options in ("-b", "--broaden"):        
        sigma = float(argument)
    elif options in ("-i", "--input"):
        infile = argument
    elif options in ("-w", "--weight"):
        weight = True 
    elif options in ("-u", "--unnorm"):
        norm = None
    else:
        assert False, "unhandled option"

if outfile==None or infile==None:
    useage()
        
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
eigenstates = numGrep("^ *([+-]?\d+\.\d+) +([+-]?\d+\.\d+e?-?\d*)",fp)
fp.close()

if norm:    
    fp = open(infile, "r")
    nKPTS = numGrep("^ *KPT *(\d+) *(.*)",fp)    
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
