#!/bin/bash

# Start a python server to drive the DFTB+ instance
sleep 2
./prerun.py &
echo "$!" > subprocess.pid
sleep 2
