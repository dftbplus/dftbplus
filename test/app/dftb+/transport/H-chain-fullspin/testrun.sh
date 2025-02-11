#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#


DFTBPLUS_CMD=$*

# Calculate the contact self energies
for contact in source drain
do
    rm -f dftb_in.hsd
    cp $contact.hsd dftb_in.hsd
    $DFTBPLUS_CMD
done

# calculate the actual device
rm -f dftb_in.hsd
cp device.hsd dftb_in.hsd
$DFTBPLUS_CMD
