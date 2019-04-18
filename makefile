# Makefile for FMM3D
#
# This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran

CC=gcc
FC=gfortran

CFLAGS=-fPIC -O3 -funroll-loops -march=native
FFLAGS = $(CFLAGS)

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS = -fopenmp
MOMPFLAGS = -lgomp -D_OPENMP

# flags for MATLAB MEX compilation..
MFLAGS=-largeArrayDims

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=mwrap

# For your OS, override the above by placing make variables in make.inc
-include make.inc

# multi-threaded libs & flags needed
ifeq ($(OMP),on)
CFLAGS += $(OMPFLAGS)
FFLAGS += $(OMPFLAGS)
endif

LIBNAME=libfmm3d
DYNAMICLIB = lib/$(LIBNAME).so
STATICLIB = lib-static/$(LIBNAME).a

# objects to compile
#
# Common objects
COM = src/Common
COMOBJS = $(COM)/besseljs3d.o $(COM)/cdjseval3d.o $(COM)/dfft.o \
	$(COM)/fmmcommon.o $(COM)/legeexps.o $(COM)/prini.o \
	$(COM)/rotgen.o $(COM)/rotproj.o $(COM)/rotviarecur.o \
	$(COM)/tree_lr_3d.o $(COM)/yrecursion.o

# Helmholtz objects
HELM = src/Helmholtz
HOBJS = $(HELM)/h3dcommon.o $(HELM)/h3dterms.o $(HELM)/h3dtrans.o \
	$(HELM)/helmrouts3d.o $(HELM)/hfmm3dpart.o $(HELM)/hfmm3dpartwrap.o \
	$(HELM)/hfmm3dpartwrap_vec.o $(HELM)/hpwrouts.o \
	$(HELM)/hwts3.o $(HELM)/numphysfour.o $(HELM)/projections.o \
	$(HELM)/quadread.o

# Laplace objects
LAP = src/Laplace
LOBJS = $(LAP)/lwtsexp_sep1.o $(LAP)/l3dterms.o $(LAP)/l3dtrans.o \
	$(LAP)/laprouts3d.o $(LAP)/lfmm3dpart.o $(LAP)/rfmm3dpartwrap.o \
	$(LAP)/rfmm3dpartwrap_vec.o $(LAP)/lwtsexp_sep2.o \
	$(LAP)/lpwrouts.o

OBJS = $(COMOBJS) $(HOBJS) $(LOBJS)

.PHONY: usage lib examples test perftest python all

default: usage

all: lib examples test perftest python

usage:
	@echo "Makefile for FMM3D. Specify what to make:"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make examples - compile and run codes in examples/"
	@echo "  make test - compile and run quick math validation tests"
	@echo "  make perftest - compile and run performance tests"
	@echo "  make python - compile and test python interfaces"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [tast] OMP=ON' for multi-threaded (otherwise single-threaded)"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
ifeq ($(OMP),ON)
	echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif
$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared $(OMPFLAGS) $(OBJS) -o $(DYNAMICLIB)


test: $(STATICLIB) test/helmrouts test/hfmm3dpart test/hfmm3dpart_vec test/laprouts test/rfmm3dpart test/rfmm3dpart_vec
	(cd test/Helmholtz; ./run_helmtest.sh)
	(cd test/Laplace; ./run_laptest.sh)

test/testhelmrouts:
	$(FC) $(FFLAGS) test/Helmholtz/test_helmrouts3d.f $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_helmrouts3d 

clean: objclean
	rm -f lib-static/*.a lib/*.so
	

objclean: 
	rm -f $(OBJS)
	rm -f test/*.o examples/*.o
