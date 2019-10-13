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