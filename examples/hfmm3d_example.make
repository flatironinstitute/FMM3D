#HOST = gcc
#HOST = gcc-openmp
#HOST = intel
HOST = intel-openmp

PROJECT = hfmm3d_example

# FC - fortran compiler
# FFLAGS - fortran compiler flags

# This make file presumes that the static library is already created
# and located at location given by 
# STATICLIB

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
    FFLAGS= -xW -O3 - xW -ip -xHost
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort
    FFLAGS= -O3 -ip -xHost -qopenmp
endif




LIBNAME=libfmm3d
STATICLIB = ../lib-static/$(LIBNAME).a

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o

.PHONY: all clean

default: all


OBJECTS = hfmm3d_example.o \
    ../src/Common/hkrand.o \
    ../src/Common/dlaran.o 

all: $(OBJECTS) 
	$(FC) $(FFLAGS)  -o $(PROJECT) $(OBJECTS) $(STATICLIB)
	./$(PROJECT)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

	

clean: 
	rm -f $(OBJECTS) $(PROJECT)
