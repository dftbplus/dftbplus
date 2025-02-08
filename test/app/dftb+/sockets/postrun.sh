#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#

############################################################################
#
#  Script for cleaning up after socket/file communication
#
############################################################################

if [ -e file.txt ]; then
    # kill anything using the file
    fuser -k $(cat file.txt)

    # Remove any left over socket file in /tmp
    rm -f $(cat file.txt)

    # remove the file itself
    rm file.txt
fi

if [ -e port.txt ]; then
    # kill anything using that port number
    fuser -n tcp -k $(cat port.txt)

    # clean up file
    rm -f port.txt
fi
