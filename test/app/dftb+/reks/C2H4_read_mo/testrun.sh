#!/bin/bash

DFTBPLUS_CMD=$*

rm -f dftb_in.hsd
cp dftb_in1.hsd dftb_in.hsd
$DFTBPLUS_CMD

rm -f dftb_in.hsd
cp dftb_in2.hsd dftb_in.hsd
$DFTBPLUS_CMD
