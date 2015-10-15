EXE=snapshot halo
default: $(EXE)
include ../Makefile.inc

snapshot: HDFLIB=
$(EXE): CXXFLAGS+=-DTEST_SIMU_IO -g -gdwarf-2 -O3
$(EXE): ../config_parser.o ../mymath.o
halo: snapshot.o
halo: CXXFLAGS+=-UTEST_SIMU_IO -DTEST_SIMU_IO_halo

