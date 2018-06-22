#!/bin/sh
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

############################################################################
#
#  Script for starting up socket/file communication
#
############################################################################

# Start a python server to drive the DFTB+ instance
sleep 2
./prerun.py &
echo "$!" > subprocess.pid
sleep 2

# run the actual calculation
$2 $1

# clean up afterwards
./postrun.sh
