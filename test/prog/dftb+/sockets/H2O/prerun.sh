#!/bin/bash

# Start a python server to drive the DFTB+ instance
./prerun.py &
echo "$!" > subprocess.pid
sleep 2
