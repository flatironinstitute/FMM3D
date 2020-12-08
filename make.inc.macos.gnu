# makefile overrides
# OS:       macOS
# Compiler: gfortran 9.X
# OpenMP:   enabled
#

CC=gcc-9
CXX=g++-9
FC=gfortran-9

ifeq ($(PREFIX),)
    FMM_INSTALL_DIR=/usr/local/lib
endif


CFLAGS += -I src 

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

# MATLAB interface:
MFLAGS += -L/usr/local/lib/gcc/9
MEX = $(shell ls -d /Applications/MATLAB_R201*.app)/bin/mex
#LIBS = -lm -lstdc++.6
#MEXLIBS= -lm -lstdc++.6 -lgfortran -ldl


