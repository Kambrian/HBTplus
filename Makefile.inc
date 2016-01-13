HDFINC=
HDFLIB=-lhdf5_hl -lhdf5

###comment out the following lines if you do not have GSL
GSLINC=
GSLINC+=-DHAS_GSL 
GSLLIB=-lgsl -lgslcblas
###comment up to here to disable GSL.

CXX=g++
CXXFLAGS+=-std=c++0x $(HDFINC) $(GSLINC) $(OMPFLAG) -O3 
CXXFLAGS+=-DENABLE_EXPERIMENTAL_PROPERTIES
<!-- CXXFLAGS+=-xHost #to get it runnable on cosma-e -->
LDLIBS+=$(HDFLIB) $(GSLLIB) $(OMPFLAG)
<!-- LINK.o=$(LINK.cpp) -->

ifeq ($(CXX), g++)
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
