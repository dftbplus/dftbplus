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
