HOST = gcc-9
HOST = gcc-9-openmp
#HOST = gcc 
#HOST = gcc-openmp

PROJECT = lfmm3d_vec_example

# CC - c compiler
# CFLAGS - fortran compiler flags

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

ifeq ($(HOST),gcc-9)
    CC=gcc-9
    CFLAGS=-std=c99 -fPIC -O3 -funroll-loops -march=native  
    CLINK= -lgfortran -lm -ldl
endif

ifeq ($(HOST),gcc-9-openmp)
    CC = gcc-9
    CFLAGS=-std=c99 -fPIC -O3 -funroll-loops -march=native -fopenmp 
    CLINK=-lgfortran -lm -ldl
endif

ifeq ($(HOST),gcc)
    CC=gcc
    CFLAGS=-std=c99 -fPIC -O3 -funroll-loops -march=native 
    CLINK= -lgfortran -lm -ldl
endif

ifeq ($(HOST),gcc-openmp)
    CC = gcc
    CFLAGS=-std=c99 -fPIC -O3 -funroll-loops -march=native -fopenmp
    CLINK= -lgfortran -lm -ldl
endif


LIBNAME=libfmm3d
STATICLIB = ../lib-static/$(LIBNAME).a

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o

.PHONY: all clean

default: all


OBJECTS = lfmm3d_vec_example.o cprini.o utils.o 

HEADERS = c/cprini.h c/utils.h c/lfmm3d_c.h     

all: $(OBJECTS) 
	$(CC) $(CFLAGS)  -o $(PROJECT) $(OBJECTS) $(CLINK) $(STATICLIB)
	./$(PROJECT)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

clean: 
	rm -f $(OBJECTS) $(PROJECT) fort.13
