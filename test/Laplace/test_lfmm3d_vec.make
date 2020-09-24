#HOST = gcc
HOST = gcc-openmp

PROJECT = int2-lfmm3d-vec

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
LAP = ../../src/Laplace

.PHONY: all clean

default: all


OBJECTS = test_lfmm3d_vec.o \
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
    $(LAP)/lfmm3d.o \
    $(LAP)/lfmm3dwrap_vec.o \
    $(LAP)/lwtsexp_sep1.o \
    $(LAP)/lwtsexp_sep2.o \
    $(LAP)/lpwrouts.o \
    $(LAP)/lndiv.o \
    $(COM)/rotproj.o \
    $(COM)/tree_routs3d.o \
    $(COM)/pts_tree3d.o \
    $(COM)/cumsum.o \
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
