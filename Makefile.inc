HDFINC=
HDFLIB=-lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5
OMPFLAG=-fopenmp

CXX=g++
CXXFLAGS+=-std=c++0x $(HDFINC) $(OMPFLAG) -g 
CXXFLAGS+=-DHBT_REAL8 -DHBT_INT8 
LDLIBS+=$(HDFLIB) $(OMPFLAG)
<!-- LINK.o=$(LINK.cpp) -->

SUBDIRS:=$(wildcard */.)
SUBDIRSCLEAN=$(addsuffix clean,$(SUBDIRS))

clean: clean_curdir $(SUBDIRSCLEAN)

clean_curdir:
	rm -rf *.o *~ core a.out $(EXE) $(OUTPUTS)

%clean: %
	$(MAKE) -C $< clean_curdir
