#!/bin/bash

# kill anything using that port number
fuser -n tcp -k $(cat port.txt)

# clean up file
rm -f port.txt
