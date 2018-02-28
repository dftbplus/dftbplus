#!/bin/bash

# kill anything using the file
fuser -k /tmp/ipi_dftb

# Remove any left over socket file in /tmp
rm -f /tmp/ipi_dftb
