#HOST = gcc
HOST = gcc-openmp

PROJECT = int2-laprouts3d

# FC - fortran compiler
# FFLAGS - fortran compiler flags

ifeq ($(HOST),gcc)
    FC=gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=x86-64 -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=x86-64 -fopenmp -std=legacy
endif

# Test objects
#
COM = ../../src/Common
LAP = ../../src/Laplace

.PHONY: all clean

default: all


OBJECTS = test_laprouts3d.o \
    $(COM)/hkrand.o \
    $(COM)/dlaran.o \
    $(COM)/prini.o \
    $(COM)/rotgen.o \
    $(COM)/legeexps.o \
    $(COM)/rotviarecur.o \
    $(COM)/yrecursion.o \
    $(LAP)/l3dterms.o \
    $(LAP)/l3dtrans.o \
    $(LAP)/laprouts3d.o \
    $(LAP)/lapkernels.o \
    $(COM)/rotproj.o \
    $(COM)/dfft.o \
    $(COM)/fmmcommon.o \

all: $(OBJECTS) 
	$(FC) $(FFLAGS)  -o $(PROJECT) $(OBJECTS) 
	./$(PROJECT)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

clean: 
	rm -f $(OBJECTS) $(PROJECT) fort.13
