OS = osx

#HOST = gcc
HOST = gcc-openmp
#HOST = intel
#HOST = intel-openmp

PROJECT = hfmm3d_example

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

ifeq ($(OS),osx)
    LDFMM = /usr/local/lib
endif

ifeq ($(OS),linux)
    LDFMM = ./../lib
endif


ifeq ($(HOST),gcc)
    FC=gfortran -L${LDFMM} 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native  
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran -L${LDFMM}
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy 
endif

ifeq ($(HOST),intel)
    FC=ifort -L${LDFMM}
    FFLAGS= -O3 -fPIC -march=native
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort -L${LDFMM}
    FFLAGS= -O3 -fPIC -march=native -qopenmp
endif

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o

.PHONY: all clean

default: all


OBJECTS = hfmm3d_vec_example.o \
    ../src/Common/hkrand.o \
    ../src/Common/dlaran.o 

all: $(OBJECTS) 
	$(FC) $(FFLAGS)  -o $(PROJECT) $(OBJECTS) -lfmm3d 
	./$(PROJECT)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

	

clean: 
	rm -f $(OBJECTS) $(PROJECT) fort.13
