
EXEC = int2

#HOST=linux-gfortran-openmp
#HOST = linux-gfortran
#HOST= linux-gfortran-openmp
#HOST = linux-ifort
HOST = macos

ifeq ($(HOST),macos)
  FC = gfortran-9
  FFLAGS = -O2 -c -w -fopenmp -march=native  
  FLINK = gfortran-9 -w -fopenmp -o $(EXEC)
  FEND = -Wl,-stack_size,0x40000000
endif

ifeq ($(HOST),linux-gfortran-openmp)
FC = gfortran
FFLAGS = -O3 -c -w --openmp 
FLINK = gfortran -w -o $(EXEC) --openmp
endif

ifeq ($(HOST),linux-gfortran)
FC = gfortran
FFLAGS = -O3 -c -w  -march=native -pg -ffast-math -ftree-vectorize -funroll-loops
FLINK = gfortran -w -o $(EXEC) -pg 
endif

ifeq ($(HOST),linux-ifort)
  FC = ifort
  FFLAGS = -O3 -c -w  -xW -ip -xHost
  FLINK = ifort -w -o $(EXEC)
endif


.PHONY: all clean list

SOURCES =  test_hfmm3d_mps.f90 \
  ../../src/Common/hkrand.f \
  ../../src/Common/dlaran.f
  # tree_lr_3d.f \
  # dlaran.f \
  # hkrand.f \
  # prini.f \
  # rotgen.f \
  # legeexps.f \
  # rotviarecur.f \
  # yrecursion.f \
  # h3dterms.f \
  # h3dtrans.f \
  # helmrouts3d.f \
  # hfmm3d.f \
  # hfmm3dwrap_vec.f \
  # hpwrouts.f \
  # hwts3.f \
  # fmmcommon.f \
  # quadread.f \
  # besseljs3d.f \
  # projections.f \
  # rotproj.f \
  # dfft.f \
  # h3dcommon.f \
  # numphysfour.f \
  # hnumphys.f \
  # hnumfour.f \
  # hwts3e.f 

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
	export DYLD_LIBRARY_PATH=../../lib:$(DYLD_LIBRARY_PATH)
	$(FLINK) $(OBJECTS) -L'../../lib-static' -lfmm3d $(FEND)
	./$(EXEC)

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



