HDFINC=
HDFLIB=-lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5

CXX=g++
CXXFLAGS+=-std=c++0x $(HDFINC) $(OMPFLAG) -O3
CXXFLAGS+=-DHBT_REAL8 -DHBT_INT8 
LDLIBS+=$(HDFLIB) $(OMPFLAG)
<!-- LINK.o=$(LINK.cpp) -->

ifeq ($(CXX), g++)
OMPFLAG=-fopenmp
else
OMPFLAG=-openmp
endif

SUBDIRS:=$(wildcard */.)
SUBDIRSCLEAN=$(addsuffix clean,$(SUBDIRS))

clean: clean_curdir $(SUBDIRSCLEAN)

clean_curdir:
	rm -rf *.o *~ core a.out $(EXE) $(OUTPUTS)

%clean: %
	$(MAKE) -C $< clean_curdir
