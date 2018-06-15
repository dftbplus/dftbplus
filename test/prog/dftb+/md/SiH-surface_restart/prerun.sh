#!/bin/bash

$1 > log.txt

tail -n 18 geo_end.xyz | cut -c 72- > velocities.dat

cp dftb_in2.hsd dftb_in.hsd
