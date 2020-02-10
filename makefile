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
CXX=g++
FC=gfortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops 



CFLAGS= -std=c99 
CFLAGS+= $(FFLAGS) 
CXXFLAGS= -std=c++11 -DSCTL_PROFILE=-1
CXXFLAGS+=$(FFLAGS)

CLIBS = -lgfortran -lm -ldl 


LIBS = -lm 

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 
MOMPFLAGS = -D_OPENMP

# flags for MATLAB MEX compilation..
MFLAGS=-largeArrayDims -DMWF77_UNDERSCORE1 
MWFLAGS=-c99complex 

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap-0.33/mwrap

MEXLIBS=-lm -lstdc++ -ldl -lgfortran

ifeq ($(FAST_KER),ON)

LIBS += -lstdc++
CLIBS += -lstdc++
OMP = ON

endif


# For your OS, override the above by placing make variables in make.inc
-include make.inc

# multi-threaded libs & flags needed
ifeq ($(OMP),ON)
CFLAGS += $(OMPFLAGS)
FFLAGS += $(OMPFLAGS)
MFLAGS += $(MOMPFLAGS)
LIBS += $(OMPLIBS)
MEXLIBS += $(OMPLIBS)
endif


LIBNAME=libfmm3d
DYNAMICLIB = $(LIBNAME).so
STATICLIB = lib-static/$(LIBNAME).a

# vectorized kernel directory
SRCDIR = ./vec-kernels/src
INCDIR = ./vec-kernels/include
LIBDIR = lib-static

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
	$(HELM)/helmrouts3d.o $(HELM)/hfmm3d.o $(HELM)/hfmm3dwrap.o \
	$(HELM)/hfmm3dwrap_legacy.o $(HELM)/hfmm3dwrap_vec.o $(HELM)/hpwrouts.o \
	$(HELM)/hwts3e.o $(HELM)/hnumphys.o $(HELM)/hnumfour.o $(HELM)/projections.o 

# Laplace objects
LAP = src/Laplace
LOBJS = $(LAP)/lwtsexp_sep1.o $(LAP)/l3dterms.o $(LAP)/l3dtrans.o \
	$(LAP)/laprouts3d.o $(LAP)/lfmm3d.o $(LAP)/lfmm3dwrap.o \
	$(LAP)/lfmm3dwrap_legacy.o $(LAP)/lfmm3dwrap_vec.o $(LAP)/lwtsexp_sep2.o \
	$(LAP)/lpwrouts.o

ifneq ($(FAST_KER),ON)
LOBJS += $(LAP)/lapkernels.o
HOBJS += $(HELM)/helmkernels.o
endif

ifeq ($(FAST_KER),ON)
LOBJS += $(LAP)/lapkernels_fast.o
HOBJS += $(HELM)/helmkernels_fast.o
COMOBJS+= $(SRCDIR)/libkernels.o
endif

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o

# C Headers and objects
COBJS = c/cprini.o c/utils.o
CHEADERS = c/cprini.h c/utils.h c/hfmm3d_c.h c/lfmm3d_c.h


OBJS = $(COMOBJS) $(HOBJS) $(LOBJS)

.PHONY: usage lib examples test perftest python all c c-examples matlab python3 big-test pw-test debug

default: usage

##all: cxxkernel lib examples test python python3 c c-examples matlab

cxxkernel: $(CXXOBJ)

$(SRCDIR)/libkernels.o: $(SRCDIR)/libkernels.cpp
		$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $^ -o $@

usage:
	@echo "Makefile for FMM3D. Specify what to make:"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make examples - compile and run fortran examples in examples/"
	@echo "  make c-examples - compile and run c examples in c/"
	@echo "  make test - compile and run validation tests (will take around 30 secs)"
	@echo "  make matlab - compile matlab interfaces"
	@echo "  make mex - generate matlab interfaces (for expert users only, requires mwrap)"
	@echo "  make python - compile and test python interfaces"
	@echo "  make python3 - compile and test python interfaces using python3"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=ON' for multi-threaded"
	@echo "  'make [task] FAST_KER=ON' for using vectorized kernel evaluation and multi-threaded (needs c++)"


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
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif
$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared -fPIC $(OMPFLAGS) $(OBJS) -o $(DYNAMICLIB) $(LIBS) 
	mv $(DYNAMICLIB) lib/


# matlab..
MWRAPFILE = fmm3d
MWRAPFILE2 = fmm3d_legacy
GATEWAY = $(MWRAPFILE)
GATEWAY2 = $(MWRAPFILE2)

matlab:	$(STATICLIB) matlab/$(GATEWAY).c matlab/$(GATEWAY2).c
	$(MEX) matlab/$(GATEWAY).c $(STATICLIB) $(MFLAGS) -output matlab/fmm3d $(MEXLIBS)
	$(MEX) matlab/$(GATEWAY2).c $(STATICLIB) $(MFLAGS) -output matlab/fmm3d_legacy $(MEXLIBS)


mex:  $(STATICLIB)
	cd matlab; $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) $(MEXLIBS); \
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY2) -mb $(MWRAPFILE2).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY2) -c $(GATEWAY2).c $(MWRAPFILE2).mw;\
	$(MEX) $(GATEWAY2).c ../$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE2) $(MEXLIBS);

#python
python: $(STATICLIB)
	cd python && export FAST_KER=$(FAST_KER) && export FLIBS='$(LIBS)' && export FFLAGS='$(FFLAGS)' && pip install -e . && cd test && pytest -s

#python
python3: $(STATICLIB)
	cd python && pip3 install -e . && cd test && python3 -m pytest -s

# testing routines
#
test: $(STATICLIB) $(TOBJS) test/helmrouts test/hfmm3d test/hfmm3d_vec test/laprouts test/lfmm3d test/lfmm3d_vec
	cd test/Helmholtz; ./test_hfmm3d
	(cd test/Helmholtz; ./run_helmtest.sh)
	(cd test/Laplace; ./run_laptest.sh)
	cat print_testreshelm.txt
	cat print_testreslap.txt
	rm print_testreshelm.txt
	rm print_testreslap.txt

test/helmrouts: 
	$(FC) $(FFLAGS) test/Helmholtz/test_helmrouts3d.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_helmrouts3d $(LIBS) 

test/hfmm3d:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_hfmm3d $(LIBS)

test/hfmm3d_vec:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_vec.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_hfmm3d_vec  $(LIBS)

test/laprouts:
	$(FC) $(FFLAGS) test/Laplace/test_laprouts3d.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/test_laprouts3d $(LIBS)

test/lfmm3d:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/test_lfmm3d $(LIBS)

test/lfmm3d_vec:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d_vec.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/test_lfmm3d_vec $(LIBS) 


test_hfmm3d_mps: $(STATICLIB) $(TOBJS)
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_mps.f90 \
  src/Helmholtz/hfmm3d_mps.f90 $(TOBJS) $(COMOBJS) $(HOBJS) \
  -o test/Helmholtz/test_hfmm3d_mps
	(cd test/Helmholtz; ./test_hfmm3d_mps)


#
##  examples
#

examples: cxxkernel $(STATICLIB) $(TOBJS) examples/ex1_helm examples/ex2_helm examples/ex3_helm \
	examples/ex1_lap examples/ex2_lap examples/ex3_lap
	./examples/lfmm3d_example
	./examples/lfmm3d_vec_example
	./examples/lfmm3d_legacy_example
	./examples/hfmm3d_example
	./examples/hfmm3d_vec_example
	./examples/hfmm3d_legacy_example

examples/ex1_lap:
	$(FC) $(FFLAGS) examples/lfmm3d_example.f $(TOBJS) $(COMOBJS) $(LOBJS) -o examples/lfmm3d_example $(LIBS)

examples/ex2_lap:
	$(FC) $(FFLAGS) examples/lfmm3d_vec_example.f $(TOBJS) $(COMOBJS) $(LOBJS) -o examples/lfmm3d_vec_example $(LIBS)

examples/ex3_lap:
	$(FC) $(FFLAGS) examples/lfmm3d_legacy_example.f $(TOBJS) $(COMOBJS) $(LOBJS) -o examples/lfmm3d_legacy_example $(LIBS)


examples/ex1_helm:
	$(FC) $(FFLAGS) examples/hfmm3d_example.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/hfmm3d_example $(LIBS)

examples/ex2_helm:
	$(FC) $(FFLAGS) examples/hfmm3d_vec_example.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/hfmm3d_vec_example $(LIBS)

examples/ex3_helm:
	$(FC) $(FFLAGS) examples/hfmm3d_legacy_example.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/hfmm3d_legacy_example $(LIBS)



# C interface
c: $(COBJS) $(OBJS) $(CHEADERS) c/lfmm3d c/hfmm3d

c/lfmm3d:
	$(CC) $(CFLAGS) c/test_lfmm3d.c $(COBJS) $(OBJS) -o c/test_lfmm3d $(CLIBS)
	time -p c/test_lfmm3d

c/hfmm3d:
	$(CC) $(CFLAGS) c/test_hfmm3d.c $(COBJS) $(OBJS) -o c/test_hfmm3d $(CLIBS)
	time -p c/test_hfmm3d



# C examples
c-examples: $(COBJS) $(OBJS) $(CHEADERS) c/ex1_lap c/ex2_lap c/ex1_helm c/ex2_helm 
	c/lfmm3d_example
	c/lfmm3d_vec_example
	c/hfmm3d_example
	c/hfmm3d_vec_example
	rm fort.13

c/ex1_lap:
	$(CC) $(CFLAGS) c/lfmm3d_example.c $(COBJS) $(OBJS) -o c/lfmm3d_example $(CLIBS)

c/ex2_lap:
	$(CC) $(CFLAGS) c/lfmm3d_vec_example.c $(COBJS) $(OBJS) -o c/lfmm3d_vec_example $(CLIBS) 

c/ex1_helm:
	$(CC) $(CFLAGS) c/hfmm3d_example.c $(COBJS) $(OBJS) -o c/hfmm3d_example $(CLIBS)

c/ex2_helm:
	$(CC) $(CFLAGS) c/hfmm3d_vec_example.c $(COBJS) $(OBJS) -o c/hfmm3d_vec_example $(CLIBS)

clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f python/*.so
	rm -rf python/build
	rm -rf python/fmm3dpy.egg-info
	rm -f examples/lfmm3d_example
	rm -f examples/lfmm3d_vec_example
	rm -f examples/lfmm3d_legacy_example
	rm -f test/Laplace/test_lfmm3d
	rm -f test/Laplace/test_lfmm3d_vec
	rm -f test/Laplace/test_lfmm3d_big
	rm -f test/Laplace/test_laprouts3d
	rm -f test/Helmholtz/test_hfmm3d
	rm -f test/Helmholtz/test_hfmm3d_vec
	rm -f test/Helmholtz/test_hfmm3d_big
	rm -f test/Helmholtz/test_helmrouts3d
	rm -f examples/hfmm3d_example
	rm -f examples/hfmm3d_vec_example
	rm -f examples/hfmm3d_legacy_example
	rm -f c/hfmm3d_example
	rm -f c/hfmm3d_vec_example
	rm -f c/lfmm3d_example
	rm -f c/lfmm3d_vec_example
	rm -f c/test_hfmm3d
	rm -f c/test_lfmm3d
	rm -f vec-kernels/src/libkernels.o

big-test: $(STATICLIB) $(TOBJS) test/test_lap_big test/test_helm_big

pw-test: $(STATICLIB) $(TOBJS) test/test_helm_pw
	./test/Helmholtz/test_hfmm3d_pw


test/test_helm_big:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_big.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_hfmm3d_big

test/test_lap_big:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d_big.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/test_lfmm3d_big

test/test_helm_pw:
	$(FC) $(FFLAGS) test/Helmholtz/test_pwrep_hfmm3d.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_hfmm3d_pw

debug: $(STATICLIB) $(TOBJS) examples/hfmm3d_deb 
	time -p examples/hfmm3d_debug


examples/hfmm3d_deb:
	$(FC) $(FFLAGS) examples/hfmm3d_debug1.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/hfmm3d_debug $(LIBS)

objclean: 
	rm -f $(OBJS) $(COBJS) $(TOBJS)
	rm -f test/*.o examples/*.o c/*.o
