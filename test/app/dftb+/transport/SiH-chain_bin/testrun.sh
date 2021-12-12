#!/bin/bash

DFTBPLUS_CMD=$*

# Calculate the contact self energies
for contact in source drain
do
    rm -f dftb_in.hsd
    cp $contact.hsd dftb_in.hsd
    $DFTBPLUS_CMD
    if [ "$contact" ==  "source" ];
    then
	echo 'AtomRange = 5 12' > Contact.$contact
    else
	echo 'AtomRange = 13 20' > Contact.$contact
    fi
    echo "Id = $contact" >> Contact.$contact
    echo "PLShiftTolerance = 1E-6" >> Contact.$contact
done

# calculate the actual device
rm -f dftb_in.hsd
cp device.hsd dftb_in.hsd
$DFTBPLUS_CMD
