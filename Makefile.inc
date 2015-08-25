CXX=icpc
CXXFLAGS+=-DHBT_REAL8 #-DHBT_INT8
LINK.o=$(LINK.cpp)

#%.o: %.cpp
#	$(CXX) -c $< $(CXXFLAGS)

SUBDIRS:=$(wildcard */.)
SUBDIRSCLEAN=$(addsuffix clean,$(SUBDIRS))

clean: clean_curdir $(SUBDIRSCLEAN)

clean_curdir:
	rm -rfv *.o *~ core a.out $(EXE)

%clean: %
	$(MAKE) -C $< clean_curdir
