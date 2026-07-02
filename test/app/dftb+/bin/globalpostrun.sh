#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

############################################################################
#
#  Script for automated testing outside of tagged files
#
############################################################################

# common text files
for FILE in output *.out *.DAT *.dat *.xyz *.gen
do
    if [ -f "$FILE" ]
    then
        if grep -iw NaN $FILE; then
            echo "Error: NaN found in file: $FILE" >&2
            echo "Error: NaN found in file: $FILE" >> stderror.log
            echo "Moving autest.tag to autotest_NaN.tag" >> stderror.log
            echo "Moving autest.tag to autotest_NaN.tag" >&2
            mv autotest.tag autotest_NaN.tag
            exit 1
        fi
    fi
done
