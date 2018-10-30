#!/bin/bash

DFTBPLUS_CMD=$*

# Calculate the contact self energies
for contact in source drain
do
    rm -f dftb_in.hsd
    cp $contact.hsd dftb_in.hsd
    $DFTBPLUS_CMD
    Ef=$(grep Fermi shiftcont_$contact.dat | sed 's/.*: *//g' | sed 's/ .*//g')
    echo -e "FermiLevel = $Ef" > Fermi.$contact
done

# calculate the actual device
rm -f dftb_in.hsd
cp device.hsd dftb_in.hsd
$DFTBPLUS_CMD
