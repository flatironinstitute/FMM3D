#!/bin/bash -e
# This is an auxiliary script to build wheels
PYTHON=$1
echo "PYTHON=$PYTHON" >> make.inc 
$PYTHON -m pip install numpy
make python-dist
