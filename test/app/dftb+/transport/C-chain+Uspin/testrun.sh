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
    if [ "$contact" ==  "source" ];
    then
	echo 'AtomRange = 9 16' > Contact.$contact
    else
	echo 'AtomRange = 17 24' > Contact.$contact
    fi
    echo "Id = $contact" >> Contact.$contact
done

# calculate the actual device
rm -f dftb_in.hsd
cp device.hsd dftb_in.hsd
$DFTBPLUS_CMD
