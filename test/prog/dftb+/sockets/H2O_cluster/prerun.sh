#!/bin/bash

# Remove any left over socket file in /tmp
rm -f /tmp/ipi_dftb_cluster

# Start a python server to drive the DFTB+ instance
./prerun.py &
echo "$!" > subprocess.pid
sleep 2
