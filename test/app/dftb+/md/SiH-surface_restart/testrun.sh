#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#


DFTBPLUS_CMD=$*

rm -f dftb_in.hsd
cp dftb_in1.hsd dftb_in.hsd
$DFTBPLUS_CMD

# process out the velocities for a restart
tail -n 18 geo_end.xyz | cut -c 72- > velocities.dat

rm -f dftb_in.hsd
cp dftb_in2.hsd dftb_in.hsd
$DFTBPLUS_CMD
