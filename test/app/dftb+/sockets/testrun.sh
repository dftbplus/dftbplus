#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

############################################################################
#
#  Script for starting up socket/file communication
#
############################################################################

DFTBPLUS_CMD=$*

# Start a python server to drive the DFTB+ instance
sleep 2
./prerun.py &
echo "$!" > subprocess.pid
sleep 2

# run the actual calculation
$DFTBPLUS_CMD

# clean up afterwards
./postrun.sh
