#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

RUN_CMD="$@"
echo "$RUN_CMD ../../../testers/test_responses"
exec $RUN_CMD ../../../testers/test_responses
