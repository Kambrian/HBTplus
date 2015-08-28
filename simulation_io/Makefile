EXE=snapshot halo
default: $(EXE)
include ../Makefile.inc

$(EXE): CXXFLAGS+=-DTEST_SIMU_IO -g -gdwarf-2
$(EXE): ../config_parser.o ../mymath.o

