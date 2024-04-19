#!/usr/bin/env bash
set -ex

DFTBPLUS_CMD=$*

./split.py

for inp in dftb_in.hsd.[0-9]*; do
  cp ${inp} dftb_in.hsd
  ${DFTBPLUS_CMD}
done
