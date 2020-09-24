#HOST = gcc
HOST = gcc-openmp

PROJECT = int2-helmrouts3d

# FC - fortran compiler
# FFLAGS - fortran compiler flags

ifeq ($(HOST),gcc)
    FC=gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy
endif

# Test objects
#
COM = ../../src/Common
HELM = ../../src/Helmholtz

.PHONY: all clean

default: all


OBJECTS = test_helmrouts3d.o \
    $(COM)/hkrand.o \
    $(COM)/dlaran.o \
    $(COM)/prini.o \
    $(COM)/rotgen.o \
    $(COM)/legeexps.o \
    $(COM)/rotviarecur.o \
    $(COM)/yrecursion.o \
    $(HELM)/h3dterms.o \
    $(HELM)/h3dtrans.o \
    $(HELM)/helmrouts3d.o \
    $(HELM)/helmkernels.o \
    $(COM)/besseljs3d.o \
    $(HELM)/projections.o \
    $(COM)/rotproj.o \
    $(COM)/dfft.o \
    $(HELM)/h3dcommon.o \
    $(COM)/fmmcommon.o \

all: $(OBJECTS) 
	$(FC) $(FFLAGS)  -o $(PROJECT) $(OBJECTS) 
	./$(PROJECT)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

clean: 
	rm -f $(OBJECTS) $(PROJECT) fort.13
