# makefile overrides
# OS:       macOS
# Compiler: gfortran X.X/Clang
# OpenMP:   enabled
#

CC=gcc
CXX=g++
FC=gfortran

FFLAGS= -fPIC -O3 -arch arm64 -std=legacy -w -mno-outline-atomics
FFLAGS_DYN= -shared -fPIC
CFLAGS= -fPIC -O3 -arch arm64 -std=c99
CXXFLAGS= -std=c++11 -DSCTL_PROFILE=-1 -fPIC -O3 -arch arm64 


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
MEX = $(shell ls -d /Applications/MATLAB_R* | sort | tail -1)/bin/mex


