#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2018  DFTB+ developers group                                  #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

############################################################################
#
#  Trivial script for moving tagged file for regression testing mechanism
#
############################################################################

DFTBPLUS_CMD=$*

# run the actual calculation
$DFTBPLUS_CMD

rm -f autotest.tag
mv results.tag autotest.tag
