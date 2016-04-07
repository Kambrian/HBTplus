HDFINC=
HDFLIB=-lhdf5_hl -lhdf5 

###comment out the following lines if you do not have GSL
GSLINC=
GSLINC+=-DHAS_GSL 
GSLLIB=-lgsl -lgslcblas
###comment up to here to disable GSL.

<!--CC=mpicc-->
<!-- CXX=g++ -->
CXX=mpiCC
CXXFLAGS+=-std=c++0x $(HDFINC) $(GSLINC) $(OMPFLAG) -O3 -Wno-deprecated -g  
CXXFLAGS+=-DENABLE_EXPERIMENTAL_PROPERTIES
<!-- CXXFLAGS+=-DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX  -idirafter /usr/include/linux #to work around intel compiler bugs -->
<!-- CXXFLAGS+=-xHost #to get it runnable on cosma-e -->
LDLIBS+=$(HDFLIB) $(OMPFLAG) $(GSLLIB) $(BOOSTLIB)
<!-- LINK.o=$(LINK.cpp) -->
LDFLAGS+=-g

ifneq (,$(filter $(CXX),mpic++ mpicxx mpiCC))
BASECXX=$(shell $(CXX) -show|cut -c1-3)
else
BASECXX=$(CXX)
endif
  
ifeq ($(BASECXX), g++)
OMPFLAG=-fopenmp
else
OMPFLAG=-openmp
endif
  
%.test.o: %.cpp
	$(CXX) $^ -DTEST_$(basename $<) $(CXXFLAGS) -c -o $@
%.test: %.test.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@
$(TEST): CXXFLAGS+=-g -gdwarf-2 -O0

SUBDIRS:=$(wildcard */.)
SUBDIRSCLEAN=$(addsuffix clean,$(SUBDIRS))

clean: clean_curdir $(SUBDIRSCLEAN)

clean_curdir:
	rm -rf *.o *~ core a.out $(EXE) $(OUTPUTS) *.test

%clean: %
	$(MAKE) -C $< clean

depend_curdir:
	makedepend --$(CXXFLAGS)-- -Y $(SRC)