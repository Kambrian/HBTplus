EXE=snapshot
default: $(EXE)
include ../Makefile.inc

snapshot: CXXFLAGS+=-DTEST_SIMU_IO -g -gdwarf-2
snapshot: ../config_parser.o ../mymath.o

