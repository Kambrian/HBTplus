EXE=simulation_io
default: $(EXE)
include ../Makefile.inc

simulation_io: CXXFLAGS+=-DTEST_SIMU_IO -g -gdwarf-2
simulation_io: ../config_parser.o ../mymath.o

