#!/usr/bin/env bash
#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2025  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#


DFTBPLUS_CMD=$*

# Calculate the contact self energies
rm -f dftb_in.hsd
cp wireBulk.hsd dftb_in.hsd
$DFTBPLUS_CMD
Ef=$(grep Fermi shiftcont_wire.dat | sed 's/.*: *//g' | sed 's/ .*//g')
echo 'AtomRange = 5 12' > Contact.wire
echo "Id = wire" >> Contact.wire
echo "FermiLevel = $Ef" >> Contact.wire
echo "PLShiftTolerance = 1E-6" >> Contact.wire

# calculate the actual device
rm -f dftb_in.hsd
cp wireEnd.hsd dftb_in.hsd
$DFTBPLUS_CMD
