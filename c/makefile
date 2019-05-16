# simple GNU makefile for FMM3D library 

# determine host platform
ifeq ($(shell uname -s), Darwin)
  HOST = osx
else
  HOST = linux
endif

COMPILERFAMILY=gcc
#COMPILERFAMILY=clang

# flags and lib locations for each platform-compiler combo...
BUILD = $(HOST)-$(COMPILERFAMILY)
ifeq ($(BUILD),linux-gcc)
  FC = gfortran
  FFLAGS = -O2 -Wall
  CC = gcc
  CFLAGS = -std=c99 -O3 -funroll-loops -march=native -fcx-limited-range
# for ubuntu...
  LINALG = -lopenblas -lblas -llapack
  LDLIBS = $(LINALG) -lm -lgfortran
  LINK = $(CC)
endif

ifeq ($(BUILD),osx-gcc)
  FC = gfortran-8
  FFLAGS = -O2 -Wall
  CC = gcc-8
  CFLAGS = -std=c99 -O3 -funroll-loops -march=native -fcx-limited-range
  LDLIBS = -framework accelerate -lm -lgfortran
  LINK = gcc-8
endif
# (could also have osx-clang,  etc ...)

.PHONY: usage lib test clean distclean
default: usage
usage:
	@echo "Makefile for fmm3d. Specify what to make:"
	@echo " make test - compile and run tests"
	@echo " make clean - remove all objects and executables made"

# location of your hellskitchen utils...
LIBOBJS = Common/cprini.o ../src/Common/d3hplratree.o ../src/Common/fmmcommon.o \
	../src/Common/dlaran.o ../src/Common/hkrand.o ../src/Common/prini.o \
	../src/Common/rotgen.o ../src/Common/rotgen2.o ../src/Common/legeexps.o \
	../src/Common/yrecursion.o ../src/Helmholtz/h3dterms.o \
	../src/Helmholtz/h3dtrans.o ../src/Helmholtz/helmrouts3d.o \
	../src/Helmholtz/hfmm3dpart.o ../src/Helmholtz/hfmm3dpartwrap.o \
	../src/Helmholtz/hpwrouts.o ../src/Helmholtz/hwts3.o \
	../src/Helmholtz/quadread.o ../src/Common/besseljs3d.o \
	../src/Helmholtz/projections.o ../src/Common/rotproj.o \
	../src/Common/dfft.o ../src/Helmholtz/h3dcommon.o \
	../src/Common/rotviarecur.o ../src/Helmholtz/numphysfour.o 
IDIRS = -ICommon -I../src

# implicit rules for objects (note -o ensures writes to same dir)
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $(IDIRS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $(IDIRS) $< -o $@
%.o: %.f90 %.h
	$(FC) -c $(FFLAGS) $(IDIRS) $< -o $@
# rule for drivers (assumed C) that have no header
%.o: %.c
	$(CC) -c $(CFLAGS) $(IDIRS) $< -o $@


# some tests and examples... (use macro to avoid retyping below)
E1 = Helmholtz/test_stos

test: $(E1) 
# make the output files appear in examples...
	($(E1)) 

$(E1): $(E1).o $(LIBOBJS)
# note libs best if come after objects:
	$(LINK) $(E1).o -o $(E1) $(LIBOBJS) $(LDLIBS)

clean:
	rm -f $(LIBOBJS) examples/int2 $(E1) 

distclean: clean
	rm -f examples/*.vtk src/fort.* examples/fort.* */*.d
