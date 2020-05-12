
EXEC = int2

#HOST=linux-gfortran-openmp
#HOST = linux-gfortran
#HOST= linux-gfortran-openmp
#HOST = linux-ifort
#HOST = macos
HOST = macos-intel-openmp

ifeq ($(HOST),macos)
  FC = gfortran-9
  FFLAGS = -O2 -c -w -fopenmp -march=native  
  FLINK = gfortran-9 -w -fopenmp -o $(EXEC)
  FEND = -Wl,-stack_size,0x40000000
endif

ifeq ($(HOST),macos-intel-openmp)
  FC = ifort
  FFLAGS = -O2 -c -w -fopenmp -march=native 
  FLINK = ifort -w -fopenmp -mkl -o $(EXEC)
  #FEND = -Wl,-stack_size,0x40000000
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -O3 -c -w --openmp 
FLINK = gfortran -w -o $(EXEC) --openmp
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -O3 -c -w -march=native -pg -ffast-math -ftree-vectorize -funroll-loops
FLINK = gfortran -w -o $(EXEC) -pg 
endif

ifeq ($(HOST),linux-ifort)
  FC = ifort
  FFLAGS = -O3 -c -w  -xW -ip -xHost
  FLINK = ifort -w -o $(EXEC)
endif

COM = ../../src/Common
HELM = ../../src/Helmholtz

.PHONY: all clean list

SOURCES =  test_hfmm3d_adjoint.f90 \
  $(COM)/hkrand.f \
  $(COM)/dlaran.f \
  $(COM)/besseljs3d.f \
  $(COM)/cdjseval3d.f \
  $(COM)/dfft.f \
  $(COM)/fmmcommon.f \
  $(COM)/legeexps.f $(COM)/prini.f \
  $(COM)/rotgen.f $(COM)/rotproj.f $(COM)/rotviarecur.f \
  $(COM)/tree_lr_3d.f $(COM)/yrecursion.f \
  $(HELM)/hfmm3d_mps.f90 $(HELM)/helmkernels.f \
  $(HELM)/h3dcommon.f $(HELM)/h3dterms.f $(HELM)/h3dtrans.f \
  $(HELM)/helmrouts3d.f $(HELM)/hfmm3d.f $(HELM)/hfmm3dwrap.f \
  $(HELM)/hfmm3dwrap_legacy.f $(HELM)/hfmm3dwrap_vec.f $(HELM)/hpwrouts.f \
  $(HELM)/hwts3e.f $(HELM)/hnumphys.f $(HELM)/hnumfour.f $(HELM)/projections.f

OBJECTS = $(patsubst %.f,%.o,$(patsubst %.f90,%.o,$(SOURCES)))

#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(EXEC)
	$(FLINK) $(OBJECTS) $(FEND)
#	export DYLD_LIBRARY_PATH=../../lib:$(DYLD_LIBRARY_PATH)
#	$(FLINK) $(OBJECTS) -L'../../lib-static' -lfmm3d $(FEND)
#	./$(EXEC)

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



