#!/bin/bash -e
# This is an auxiliary script to test wheels
PYTHON=$1
$PYTHON -m pip install pytest
$PYTHON -m pip install fmm3dpy -f wheelhouse/ --no-index --no-cache
$PYTHON -m pytest -s python/test
