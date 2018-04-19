#!/bin/bash

# kill anything using the file
fuser -k $(cat file.txt)

# Remove any left over socket file in /tmp
rm -f $(cat file.txt)

# remove the file itself
rm file.txt
