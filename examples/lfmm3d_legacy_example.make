#HOST = gcc
HOST = gcc-openmp
#HOST = intel
#HOST = intel-openmp

PROJECT = lfmm3d_vec_example

# FC - fortran compiler
# FFLAGS - fortran compiler flags

# This make file presumes that the static library is already created
# and located at location given by 
# STATICLIB
#
# It also assumes that the static library is compiled using
# the same compiler that you are using to run the make file with
#
# In case you wish to you use a different compiler, make sure to include 
# the cross compiled libraries. See fmm3d.readthedocs.io/install.html
# for additional info

ifeq ($(HOST),gcc)
    FC=gfortran
    FFLAGS=-fPIC -O3 -funroll-loops -march=native  
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp 
endif

ifeq ($(HOST),intel)
    FC=ifort
    FFLAGS= -O3 -xW -ip -xHost
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort
    FFLAGS= -O3 -xW -ip -xHost -qopenmp
endif




LIBNAME=libfmm3d
STATICLIB = ../lib-static/$(LIBNAME).a

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o

.PHONY: all clean

default: all


OBJECTS = lfmm3d_legacy_example.o \
    ../src/Common/hkrand.o \
    ../src/Common/dlaran.o 

all: $(OBJECTS) 
	$(FC) $(FFLAGS)  -o $(PROJECT) $(OBJECTS) $(STATICLIB)
	./$(PROJECT)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

	

clean: 
	rm -f $(OBJECTS) $(PROJECT) fort.13
