#!/bin/bash

grep -IUlr --exclude-dir=.git $'\r' .
