#!/usr/bin/env bash

grep -IUlr --exclude-dir=.git $'\r' .
