# makefile overrides
# OS:       macOS
# Compiler: gfortran X.X
# OpenMP:   enabled
#

CC=gcc
CXX=g++
FC=gfortran

ifeq ($(PREFIX),)
    FMM_INSTALL_DIR=/usr/local/lib
endif


CFLAGS += -I src 

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

# MATLAB interface:
FDIR=$$(dirname `gfortran --print-file-name libgfortran.dylib`)
MFLAGS +=-L${FDIR}
MEX = $(shell ls -d /Applications/MATLAB_R20**.app)/bin/mex
#LIBS = -lm -lstdc++.6
#MEXLIBS= -lm -lstdc++.6 -lgfortran -ldl


