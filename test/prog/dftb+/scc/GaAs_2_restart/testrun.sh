#!/bin/bash

rm -f dftb_in.hsd
cp dftb_in1.hsd dftb_in.hsd
$1 > log.txt

rm -f dftb_in.hsd
cp dftb_in2.hsd dftb_in.hsd
$2 $1
