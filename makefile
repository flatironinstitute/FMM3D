# Makefile for FMM3D
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 


# compiler, and linking from C, fortran
CC=gcc
CXX=g++
FC=gfortran


# set compiler flags for c and fortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy 
CFLAGS= -fPIC -O3 -march=native -funroll-loops -std=c99 
CXXFLAGS= -std=c++11 -DSCTL_PROFILE=-1 -fPIC -O3 -march=native -funroll-loops

# set linking libraries
CLIBS = -lgfortran -lm -ldl 
LIBS = -lm

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

# Python Exetucable
PYTHON=python


# flags for MATLAB MEX compilation..
MFLAGS=-compatibleArrayDims -DMWF77_UNDERSCORE1 
MWFLAGS=-c99complex 
MOMPFLAGS = -D_OPENMP

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran

FMM_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	FMM_INSTALL_DIR = ${HOME}/lib
endif

DYLIBS = $(LIBS)

LIBNAME=libfmm3d
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LLINKLIB = -lfmm3d



# For your OS, override the above by placing make variables in make.inc
-include make.inc

# additional compile flags for FAST_KER
ifeq ($(FAST_KER),ON)
  LIBS += -lstdc++
  DYLIBS += -lstdc++
  CLIBS += -lstdc++
  FFLAGS += -lstdc++
  CFLAGS += -lstdc++
  OMP = ON
endif


# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)
  CFLAGS += $(OMPFLAGS)
  FFLAGS += $(OMPFLAGS)
  MFLAGS += $(MOMPFLAGS)
  LIBS += $(OMPLIBS)
  DYLIBS += $(OMPLIBS)
  MEXLIBS += $(OMPLIBS)
endif


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
	$(COM)/tree_routs3d.o $(COM)/pts_tree3d.o $(COM)/yrecursion.o \
	$(COM)/cumsum.o

# Helmholtz objects
HELM = src/Helmholtz
HOBJS = $(HELM)/h3dcommon.o $(HELM)/h3dterms.o $(HELM)/h3dtrans.o \
	$(HELM)/helmrouts3d.o $(HELM)/hfmm3d.o $(HELM)/hfmm3dwrap.o \
	$(HELM)/hfmm3dwrap_legacy.o $(HELM)/hfmm3dwrap_vec.o $(HELM)/hpwrouts.o \
	$(HELM)/hwts3e.o $(HELM)/hnumphys.o $(HELM)/hnumfour.o $(HELM)/projections.o \
	$(HELM)/hfmm3d_mps.o $(HELM)/hfmm3d_memest.o

# Laplace objects
LAP = src/Laplace
LOBJS = $(LAP)/lwtsexp_sep1.o $(LAP)/l3dterms.o $(LAP)/l3dtrans.o \
	$(LAP)/laprouts3d.o $(LAP)/lfmm3d.o $(LAP)/lfmm3dwrap.o \
	$(LAP)/lfmm3dwrap_legacy.o $(LAP)/lfmm3dwrap_vec.o $(LAP)/lwtsexp_sep2.o \
	$(LAP)/lpwrouts.o

# Stokes objects
STOK = src/Stokes
STOBJS = $(STOK)/stfmm3d.o $(STOK)/stokkernels.o

# Maxwell objects
EM = src/Maxwell
EMOBJS = $(EM)/emfmm3d.o

ifneq ($(FAST_KER),ON)
LOBJS += $(LAP)/lapkernels.o
LOBJS += $(LAP)/lndiv.o
HOBJS += $(HELM)/helmkernels.o
HOBJS += $(HELM)/hndiv.o
endif

ifeq ($(FAST_KER),ON)
LOBJS += $(LAP)/lapkernels_fast.o
LOBJS += $(LAP)/lndiv_fast.o
HOBJS += $(HELM)/helmkernels_fast.o
HOBJS += $(HELM)/hndiv_fast.o
COMOBJS+= $(SRCDIR)/libkernels.o
endif

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o

# C Headers and objects
COBJS = c/cprini.o c/utils.o
CHEADERS = c/cprini.h c/utils.h c/hfmm3d_c.h c/lfmm3d_c.h


OBJS = $(COMOBJS) $(HOBJS) $(LOBJS) $(STOBJS) $(EMOBJS)

.PHONY: usage install lib examples test test-ext python all c c-examples matlab big-test pw-test debug test-dyn matlab-dyn python-dyn python-dist

default: usage

cxxkernel: $(CXXOBJ)

$(SRCDIR)/libkernels.o: $(SRCDIR)/libkernels.cpp
		$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $^ -o $@

usage:
	@echo "Makefile for FMM3D. Specify what to make:"
	@echo "  make install - compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR) - compile and install the main library at custom location given by PREFIX"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test - compile and run validation tests (will take a couple of mins)"
	@echo "  make matlab - compile matlab interfaces"
	@echo "  make python - compile and test python interfaces"
	@echo "  make test-dyn - test successful installation by validation tests linked to dynamic library (will take a couple of mins)"
	@echo "  make matlab-dyn - compile matlab interfaces with dynamic library linking"
	@echo "  make python-dyn - compile and test python interfaces with dynamic library linking"
	@echo "  make python-dist - compile python interfaces for distribution"
	@echo "  make examples - compile and run fortran examples in examples/"
	@echo "  make c-examples - compile and run c examples in c/"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo "  make mex - generate matlab interfaces (for expert users only, requires mwrap)"
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=OFF' for single-threaded"
	@echo "  'make [task] FAST_KER=ON' for using vectorized kernel evaluation and multi-threaded (needs c++)"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif


install: $(STATICLIB) $(DYNAMICLIB)
	echo $(FMM_INSTALL_DIR)
	mkdir -p $(FMM_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(FMM_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(FMM_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp lib/$(LIMPLIB) $(FMM_INSTALL_DIR)/
	@echo "Make sure to include " $(FMM_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(FMM_INSTALL_DIR) " -lfmm3d"


$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/
$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared -fPIC $(OBJS) -o $(DYNAMICLIB) $(DYLIBS) 
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/


# matlab..
MWRAPFILE = fmm3d
MWRAPFILE2 = fmm3d_legacy
GATEWAY = $(MWRAPFILE)
GATEWAY2 = $(MWRAPFILE2)

matlab:	$(STATICLIB) matlab/$(GATEWAY).c matlab/$(GATEWAY2).c
	$(MEX) matlab/$(GATEWAY).c lib-static/$(STATICLIB) $(MFLAGS) \
	-output matlab/fmm3d $(MEXLIBS) 
	$(MEX) matlab/$(GATEWAY2).c lib-static/$(STATICLIB) $(MFLAGS) \
	-output matlab/fmm3d_legacy $(MEXLIBS) 

matlab-dyn:	$(DYNAMICLIB) matlab/$(GATEWAY).c matlab/$(GATEWAY2).c
	$(MEX) matlab/$(GATEWAY).c $(MFLAGS) \
	-output matlab/fmm3d $(MEXLIBS) -L$(FMM_INSTALL_DIR) $(LLINKLIB) 
	$(MEX) matlab/$(GATEWAY2).c $(MFLAGS) \
	-output matlab/fmm3d_legacy $(MEXLIBS) -L$(FMM_INSTALL_DIR) $(LLINKLIB)

mex:  $(STATICLIB)
	cd matlab; $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) \
	$(MEXLIBS); \
	$(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY2) -mb \
	$(MWRAPFILE2).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY2) -c $(GATEWAY2).c $(MWRAPFILE2).mw;\
	$(MEX) $(GATEWAY2).c ../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE2) \
	$(MEXLIBS);

#python
python: $(STATICLIB)
	cd python && \
	FMM_FLIBS='$(LIBS) $(OMPFLAGS)' $(PYTHON) -m pip install -e . && \
	$(PYTHON) -m pytest test/ -s

python-dist: $(STATICLIB)
	cd python && \
	FMM_FLIBS='$(LIBS) $(OMPFLAGS)' $(PYTHON) setup.py bdist_wheel


# testing routines
#
test: $(STATICLIB) $(TOBJS) test/helmrouts test/hfmm3d test/hfmm3d_vec test/hfmm3d_scale test/laprouts test/lfmm3d test/lfmm3d_vec test_hfmm3d_mps test/lfmm3d_scale test/stfmm3d test/stokkernels test/emfmm3d 
	(cd test/Helmholtz; ./run_helmtest.sh)
	(cd test/Laplace; ./run_laptest.sh)
	(cd test/Stokes; ./run_stoktest.sh)
	(cd test/Maxwell; ./run_emtest.sh)
	cat print_testreshelm.txt
	cat print_testreslap.txt
	cat print_testresstok.txt
	cat print_testresem.txt
	rm print_testreshelm.txt
	rm print_testreslap.txt
	rm print_testresstok.txt
	rm print_testresem.txt

test-dyn: $(DYNAMICLIB) $(TOBJS) test/helmrouts-dyn test/hfmm3d-dyn test/hfmm3d_vec-dyn test/hfmm3d_scale-dyn test/laprouts-dyn test/lfmm3d-dyn test/lfmm3d_vec-dyn test_hfmm3d_mps-dyn test/lfmm3d_scale-dyn
	(cd test/Helmholtz; ./run_helmtest.sh)
	(cd test/Laplace; ./run_laptest.sh)
	cat print_testreshelm.txt
	cat print_testreslap.txt
	rm print_testreshelm.txt
	rm print_testreslap.txt

test-ext: $(STATICLIB) $(TOBJS) test/helmrouts test/hfmm3d test/hfmm3d_vec test/hfmm3d_zkbig test/hfmm3d_scale test/laprouts test/lfmm3d test/lfmm3d_vec test_hfmm3d_mps test/lfmm3d_scale
	(cd test/Helmholtz; ./run_helmtest_ext.sh)
	(cd test/Laplace; ./run_laptest.sh)
	cat print_testreshelm.txt
	cat print_testreslap.txt
	rm print_testreshelm.txt
	rm print_testreslap.txt

test/helmrouts: 
	$(FC) $(FFLAGS) test/Helmholtz/test_helmrouts3d.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/int2-test-helmrouts3d $(LIBS) 

test/hfmm3d:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/int2-test-hfmm3d $(LIBS)

test/hfmm3d_zkbig:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_zkbig.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/int2-test-hfmm3d-zkbig $(LIBS)

test/hfmm3d_scale:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_scale.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/int2-test-hfmm3d-scale $(LIBS)

test/hfmm3d_vec:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_vec.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/int2-test-hfmm3d-vec  $(LIBS)

test/laprouts:
	$(FC) $(FFLAGS) test/Laplace/test_laprouts3d.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/int2-test-laprouts3d $(LIBS)

test/lfmm3d:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/int2-test-lfmm3d $(LIBS)

test/lfmm3d_scale:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d_scale.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/int2-test-lfmm3d-scale $(LIBS)

test/lfmm3d_vec:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d_vec.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/int2-test-lfmm3d-vec $(LIBS) 

test/stfmm3d: 
	$(FC) $(FFLAGS) test/Stokes/test_stfmm3d.f $(TOBJS) $(COMOBJS) $(LOBJS) $(STOBJS) -o test/Stokes/int2-test-stfmm3d $(LIBS) 

test/stokkernels: 
	$(FC) $(FFLAGS) test/Stokes/test_stokkernels.f $(TOBJS) $(COMOBJS) $(LOBJS) $(STOBJS) -o test/Stokes/int2-test-stokkernels $(LIBS)

test/emfmm3d: 
	$(FC) $(FFLAGS) test/Maxwell/test_emfmm3d.f $(TOBJS) $(COMOBJS) $(HOBJS) $(EMOBJS) -o test/Maxwell/int2-test-emfmm3d $(LIBS) 

test_hfmm3d_mps: $(STATICLIB) $(TOBJS)
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_mps.f90 \
  $(TOBJS) $(COMOBJS) $(HOBJS) $(LIBS)\
  -o test/Helmholtz/int2-test-hfmm3d-mps


## Linking against dynamic libraries
#
#
test/helmrouts-dyn: 
	$(FC) $(FFLAGS) test/Helmholtz/test_helmrouts3d.f $(TOBJS) -o test/Helmholtz/int2-test-helmrouts3d -L$(FMM_INSTALL_DIR) $(LLINKLIB) 

test/hfmm3d-dyn:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d.f $(TOBJS) -o test/Helmholtz/int2-test-hfmm3d -L$(FMM_INSTALL_DIR) $(LLINKLIB)  

test/hfmm3d_scale-dyn:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_scale.f $(TOBJS) -o test/Helmholtz/int2-test-hfmm3d-scale -L$(FMM_INSTALL_DIR) $(LLINKLIB) 

test/hfmm3d_vec-dyn:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_vec.f $(TOBJS) -o test/Helmholtz/int2-test-hfmm3d-vec -L$(FMM_INSTALL_DIR) $(LLINKLIB) 

test/laprouts-dyn:
	$(FC) $(FFLAGS) test/Laplace/test_laprouts3d.f $(TOBJS) -o test/Laplace/int2-test-laprouts3d -L$(FMM_INSTALL_DIR) $(LLINKLIB) 

test/lfmm3d-dyn:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d.f $(TOBJS) -o test/Laplace/int2-test-lfmm3d -L$(FMM_INSTALL_DIR) $(LLINKLIB) 

test/lfmm3d_scale-dyn:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d_scale.f $(TOBJS) -o test/Laplace/int2-test-lfmm3d-scale -L$(FMM_INSTALL_DIR) $(LLINKLIB)
	 
test/lfmm3d_vec-dyn:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d_vec.f $(TOBJS) -o test/Laplace/int2-test-lfmm3d-vec -L$(FMM_INSTALL_DIR) $(LLINKLIB) 


test_hfmm3d_mps-dyn: $(DYNAMICLIB) $(TOBJS)
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_mps.f90 \
  $(TOBJS) \
  -o test/Helmholtz/int2-test-hfmm3d-mps -L$(FMM_INSTALL_DIR) $(LLINKLIB)


#
##  examples
#

examples: cxxkernel $(STATICLIB) $(TOBJS) examples/ex1_helm examples/ex2_helm examples/ex3_helm \
	examples/ex1_lap examples/ex2_lap examples/ex3_lap
	./examples/int2-lfmm3d-example
	./examples/int2-lfmm3d-vec-example
	./examples/int2-lfmm3d-legacy-example
	./examples/int2-hfmm3d-example
	./examples/int2-hfmm3d-vec-example
	./examples/int2-hfmm3d-legacy-example

examples/ex1_lap:
	$(FC) $(FFLAGS) examples/lfmm3d_example.f $(TOBJS) $(COMOBJS) $(LOBJS) -o examples/int2-lfmm3d-example $(LIBS)

examples/ex2_lap:
	$(FC) $(FFLAGS) examples/lfmm3d_vec_example.f $(TOBJS) $(COMOBJS) $(LOBJS) -o examples/int2-lfmm3d-vec-example $(LIBS)

examples/ex3_lap:
	$(FC) $(FFLAGS) examples/lfmm3d_legacy_example.f $(TOBJS) $(COMOBJS) $(LOBJS) -o examples/int2-lfmm3d-legacy-example $(LIBS)


examples/ex1_helm:
	$(FC) $(FFLAGS) examples/hfmm3d_example.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/int2-hfmm3d-example $(LIBS)

examples/ex2_helm:
	$(FC) $(FFLAGS) examples/hfmm3d_vec_example.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/int2-hfmm3d-vec-example $(LIBS)

examples/ex3_helm:
	$(FC) $(FFLAGS) examples/hfmm3d_legacy_example.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/int2-hfmm3d-legacy-example $(LIBS)



# C interface
c: $(COBJS) $(OBJS) $(CHEADERS) c/lfmm3d c/hfmm3d

c/lfmm3d:
	$(CC) $(CFLAGS) c/test_lfmm3d.c $(COBJS) $(OBJS) -o c/int2-test-lfmm3d $(CLIBS)
	time -p c/int2-test-lfmm3d

c/hfmm3d:
	$(CC) $(CFLAGS) c/test_hfmm3d.c $(COBJS) $(OBJS) -o c/int2-test-hfmm3d $(CLIBS)
	time -p c/int2-test-hfmm3d



# C examples
c-examples: $(COBJS) $(OBJS) $(CHEADERS) c/ex1_lap c/ex2_lap c/ex1_helm c/ex2_helm 
	c/int2-lfmm3d-example
	c/int2-lfmm3d-vec-example
	c/int2-hfmm3d-example
	c/int2-hfmm3d-vec-example
	rm fort.13

c/ex1_lap:
	$(CC) $(CFLAGS) c/lfmm3d_example.c $(COBJS) $(OBJS) -o c/int2-lfmm3d-example $(CLIBS)

c/ex2_lap:
	$(CC) $(CFLAGS) c/lfmm3d_vec_example.c $(COBJS) $(OBJS) -o c/int2-lfmm3d-vec-example $(CLIBS) 

c/ex1_helm:
	$(CC) $(CFLAGS) c/hfmm3d_example.c $(COBJS) $(OBJS) -o c/int2-hfmm3d-example $(CLIBS)

c/ex2_helm:
	$(CC) $(CFLAGS) c/hfmm3d_vec_example.c $(COBJS) $(OBJS) -o c/int2-hfmm3d-vec-example $(CLIBS)

clean: objclean
	rm -f lib-static/*.a lib/*.so lib/*.dll lib/*.lib
	rm -f python/fmm3dpy*.so
	rm -rf python/build
	rm -rf python/dist
	rm -rf python/fmm3dpy.egg-info
	rm -f examples/lfmm3d_example
	rm -f examples/lfmm3d_vec_example
	rm -f examples/lfmm3d_legacy_example
	rm -f test/Laplace/int2-*
	rm -f test/Helmholtz/int2-*
	rm -f examples/hfmm3d_example
	rm -f examples/hfmm3d_vec_example
	rm -f examples/hfmm3d_legacy_example
	rm -f c/int2-*
	rm -f vec-kernels/src/libkernels.o

big-test: $(STATICLIB) $(TOBJS) test/test_lap_big test/test_helm_big

pw-test: $(STATICLIB) $(TOBJS) test/test_helm_pw
	./test/Helmholtz/test_hfmm3d_pw

test/test_helm_big:
	$(FC) $(FFLAGS) test/Helmholtz/test_hfmm3d_big.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_hfmm3d_big $(LIBS)

test/test_lap_big:
	$(FC) $(FFLAGS) test/Laplace/test_lfmm3d_big.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/Laplace/test_lfmm3d_big $(LIBS)

test/test_helm_pw:
	$(FC) $(FFLAGS) test/Helmholtz/test_pwrep_hfmm3d.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/Helmholtz/test_hfmm3d_pw

debug: $(STATICLIB) $(TOBJS) examples/lfmm3d_deb 
	time -p examples/lfmm3d_debug

examples/hfmm3d_deb:
	$(FC) $(FFLAGS) examples/hfmm3d_debug1.f $(TOBJS) $(COMOBJS) $(HOBJS) -o examples/hfmm3d_debug $(LIBS)

examples/lfmm3d_deb:
	$(FC) $(FFLAGS) examples/lfmm3d_debug1.f $(TOBJS) $(COMOBJS) $(LOBJS) -o examples/lfmm3d_debug $(LIBS)



objclean: 
	rm -f $(OBJS) $(COBJS) $(TOBJS)
	rm -f test/*.o examples/*.o c/*.o
