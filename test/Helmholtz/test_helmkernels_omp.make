#HOST = gcc
HOST = gcc-openmp
#HOST = nag

PROJECT = int2-helmkernels-omp

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

ifeq ($(HOST),nag)
    FC=nagfor
    FFLAGS=-PIC -O3 -Ounroll=1 -f90_sign -dcfuns -dusty -w=x77 -w=unreffed -w=unused -ieee=full -openmp
    OMPFLAGS= -openmp
    OMPLIBS = -lf71omp64 
    CLIBS = -lm -ldl
endif

# Test objects
#
COM = ../../src/Common
HELM = ../../src/Helmholtz

.PHONY: all clean

default: all


OBJECTS = test_helmkernels_omp.o ../../src/Helmholtz/helmkernels.o

all: $(OBJECTS) 
	$(FC) $(FFLAGS)  -o $(PROJECT) $(OBJECTS) 
	./$(PROJECT)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

clean: 
	rm -f $(OBJECTS) $(PROJECT) fort.13
