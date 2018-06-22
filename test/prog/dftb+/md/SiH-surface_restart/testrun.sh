#!/bin/bash

rm -f dftb_in.hsd
cp dftb_in1.hsd dftb_in.hsd
$1 > log.txt

# process out the velocities for a restart
tail -n 18 geo_end.xyz | cut -c 72- > velocities.dat

rm -f dftb_in.hsd
cp dftb_in2.hsd dftb_in.hsd
$2 $1
