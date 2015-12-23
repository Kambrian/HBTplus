HDFINC=
HDFLIB=-lhdf5_hl -lhdf5 

<!--CC=mpicc-->
<!-- CXX=g++ -->
CXX=mpicxx
CXXFLAGS+=-std=c++0x $(HDFINC) $(OMPFLAG) -O3 -Wno-deprecated -g  
<!-- CXXFLAGS+=-DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX  -idirafter /usr/include/linux #to work around intel compiler bugs -->
<!-- CXXFLAGS+=-xHost #to get it runnable on cosma-e -->
LDLIBS+=$(HDFLIB) $(OMPFLAG) $(BOOSTLIB)
<!-- LINK.o=$(LINK.cpp) -->

ifeq ($(CXX), mpic++)
BASECXX=$(shell mpic++ -show|cut -c1-3)
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
