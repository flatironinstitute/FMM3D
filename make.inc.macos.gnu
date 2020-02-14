# makefile overrides
# OS:       macOS
# Compiler: gfortran 9.X
# OpenMP:   enabled
#

CC=gcc-9
CXX=g++-9
FC=gfortran-9
FFLAGS= -fPIC -O3 -march=native -funroll-loops -lstdc++

CFLAGS += -I src 
CLINK += -Wl,-stack_size,0x40000000

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

# MATLAB interface:
MFLAGS += -L/usr/local/lib/gcc/9
MEX = $(shell ls -d /Applications/MATLAB_R201*.app)/bin/mex
#LIBS = -lm -lstdc++.6
#MEXLIBS= -lm -lstdc++.6 -lgfortran -ldl


