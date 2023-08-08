AR = ar cru

FC = gfortran-9
#FC = ifort
FFLAGS = -Ofast -fopenmp -march=native -Wall -unroll-aggressive
#FFLAGS = -Ofast -unroll-aggressive -xCORE-AVX512 -qopt-report:5 -qopt-zmm-usage=high -fopenmp

CXX = g++-9
#CXX = icpc
CXXFLAGS = -Ofast -fopenmp -march=native -Wall -std=c++11 -unroll-aggressive # need C++11
#CXXFLAGS = -Ofast -unroll-aggressive -xCORE-AVX512 -qopt-report:5 -qopt-zmm-usage=high -fopenmp -Wall -std=c++11# need C++11

#CXXFLAGS += -DSCTL_HAVE_SVML -mkl  # enable SVML with FC=ifort and CXX=icpc

#CXXFLAGS += -DSCTL_MEMDEBUG # Enable memory checks
CXXFLAGS += -DSCTL_PROFILE=5 -DSCTL_VERBOSE # Enable profiling




LDLIBS = -L$(LIBDIR) -lkernels -lstdc++

RM = rm -f
MKDIRS = mkdir -p

INCDIR = ./include
SRCDIR = ./src
OBJDIR = ./obj
LIBDIR = ./lib
BINDIR = ./bin

TARGET_LIB = $(LIBDIR)/libkernels.a

TARGET_BIN = \
						 $(BINDIR)/test-fortran \
						 $(BINDIR)/test-cpp

all : $(TARGET_LIB) $(TARGET_BIN)

$(BINDIR)/%: $(SRCDIR)/%.f90 $(TARGET_LIB)
	-@$(MKDIRS) $(dir $@)
	$(FC) $(FFLAGS) -I$(INCDIR) $^ $(LDLIBS) -o $@

$(BINDIR)/%: $(SRCDIR)/%.cpp $(TARGET_LIB)
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) $^ $(LDLIBS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $^ -o $@

$(LIBDIR)/%.a: $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(AR) $@ $^
	ranlib $@

clean:
	$(RM) -r $(BINDIR)/* $(LIBDIR)/* $(OBJDIR)/*

