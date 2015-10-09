EXE=tree
default: $(EXE)
include ../Makefile.inc

$(EXE): CXXFLAGS+=-DTEST_TREE -g -gdwarf-2
$(EXE): ../config_parser.o ../mymath.o ../simulation_io/snapshot.o ../simulation_io/halo.o


