CXX=g++
CXXFLAGS+=-std=c++0x
CXXFLAGS+=-DHBT_REAL8 -DHBT_INT8
LINK.o=$(LINK.cpp)

SUBDIRS:=$(wildcard */.)
SUBDIRSCLEAN=$(addsuffix clean,$(SUBDIRS))

clean: clean_curdir $(SUBDIRSCLEAN)

clean_curdir:
	rm -rf *.o *~ core a.out $(EXE)

%clean: %
	$(MAKE) -C $< clean_curdir
