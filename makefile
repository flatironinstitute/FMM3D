# Makefile for FMM3D.

# This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it hard to
# stay up to date with the repo version). Rather, in order to change
# OS/environment-specific compilers and flags, create the file make.inc, which
# overrides the defaults below (which are for ubuntu linux/GCC system).
# See docs/install.rst, and make.inc.*

# compilers, and linking from C, fortran...
CXX=g++
CC=gcc
FC=gfortran

# compile flags for GCC, baseline single-threaded, double precision case...
# Notes: 1) -Ofast breaks isfinite() & isnan(), so use -O3 which now is as fast
#        2) -fcx-limited-range for fortran-speed complex arith in C++
CFLAGS   = -fPIC -O3 -funroll-loops -march=native -fcx-limited-range
# tell examples where to find header files...
CFLAGS   += -I src
FFLAGS   = $(CFLAGS)

# decide name of obj files and finufft library we're building...
LIBNAME=lib_hfmm3d
DYNAMICLIB = lib/$(LIBNAME).so
STATICLIB = lib-static/$(LIBNAME).a

# ======================================================================

# objects to compile: common...
COBJS = src/Common/d3hplratree.o src/Common/dfft.o src/Common/fmmcommon.o src/Common/legeexps.o \
	src/Common/rotgen.o src/Common/rotproj.o src/Common/rotviarecur.o src/Common/yrecursion.o
# just the different kernels separately...
HOBJS = $(COBJS) src/Helmholtz/h3dterms.o src/Helmholtz/helmrouts3d.o src/Helmholtz/hfmm3dmain.o \
	src/Helmholtz/hfmm3dpart.o src/Helmholtz/hpwdrouts.o src/Helmholtz/hwts3.o \
	src/Helmholtz/projections.o 

HEADERS = src/spreadinterp.h src/finufft.h src/dirft.h src/common.h src/defs.h src/utils.h fortran/finufft_f.h

.PHONY: usage lib test 

default: usage

all: lib test

usage:
	@echo "Makefile for FMM3D library. Specify what to make:"
	@echo " make lib - compile the main library (in lib/ and lib-static/)"
	@echo " make test - compile and run quick math validation tests"
	@echo " make objclean - remove all object files, preserving lib & MEX"
	@echo "For faster (multicore) making, append the flag -j"
	@echo ""
	@echo ""
	@echo "Also see docs/install.rst"

# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@

# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
ifeq ($(OMP),OFF)
	echo "$(STATICLIB) and $(DYNAMICLIB) built, single-thread versions"
else
	echo "$(STATICLIB) and $(DYNAMICLIB) built, multithreaded versions"
endif
$(STATICLIB): $(OBJS) $(HEADERS)
	ar rcs $(STATICLIB) $(OBJS)
$(DYNAMICLIB): $(OBJS) $(HEADERS)
	$(CXX) -shared $(OMPFLAGS) $(LIBSFFT) $(OBJS) -o $(DYNAMICLIB)
# here $(OMPFLAGS) $(LIBSFFT) is needed for mac osx.
# see: http://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html

# validation tests... (most link to .o allowing testing pieces separately)
test: $(STATICLIB) test/testutils test/finufft1d_test test/finufft2d_test test/finufft3d_test test/dumbinputs test/finufft2dmany_test
	(cd test; \
	export FINUFFT_REQ_TOL=$(REQ_TOL); \
	export FINUFFT_CHECK_TOL=$(CHECK_TOL); \
	./check_finufft.sh)
test/testutils: test/testutils.cpp src/utils.o src/utils.h $(HEADERS)
	$(CXX) $(CXXFLAGS) test/testutils.cpp src/utils.o -o test/testutils
test/finufft1d_test: test/finufft1d_test.cpp $(OBJS1) $(HEADERS)
	$(CXX) $(CXXFLAGS) test/finufft1d_test.cpp $(OBJS1) $(LIBSFFT) -o test/finufft1d_test
test/finufft2d_test: test/finufft2d_test.cpp $(OBJS2) $(HEADERS)
	$(CXX) $(CXXFLAGS) test/finufft2d_test.cpp $(OBJS2) $(LIBSFFT) -o test/finufft2d_test
test/finufft3d_test: test/finufft3d_test.cpp $(OBJS3) $(HEADERS)
	$(CXX) $(CXXFLAGS) test/finufft3d_test.cpp $(OBJS3) $(LIBSFFT) -o test/finufft3d_test
test/dumbinputs: test/dumbinputs.cpp $(STATICLIB) $(HEADERS)
	$(CXX) $(CXXFLAGS) test/dumbinputs.cpp $(STATICLIB) $(LIBSFFT) -o test/dumbinputs
test/finufft2dmany_test: test/finufft2dmany_test.cpp $(OBJS2) $(HEADERS)
	$(CXX) $(CXXFLAGS) test/finufft2dmany_test.cpp $(OBJS2) $(LIBSFFT) -o test/finufft2dmany_test

