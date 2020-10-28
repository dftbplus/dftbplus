#!/bin/bash

DFTBPLUS_CMD=$*

for i in step1 step2
do
    cp $i.hsd dftb_in.hsd
    $DFTBPLUS_CMD
done
