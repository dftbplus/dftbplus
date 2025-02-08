#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

set -ex

DFTBPLUS_CMD=$*

./split.py

for inp in dftb_in.hsd.[0-9]*; do
  cp ${inp} dftb_in.hsd
  ${DFTBPLUS_CMD}
done
