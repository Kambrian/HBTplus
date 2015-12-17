SRC=$(wildcard *.cpp)
EXE=$(basename $(SRC))
default: $(EXE)
test_parse:  ../config_parser.o #../mymath.o
test_closestfactor: ../mymath.o

OUTPUTS=$(wildcard *.hdf5)

include ../../Makefile.inc
test_omptask:CXXFLAGS+=-O0
test_boostmpi_mem: CXXFLAGS+=-DUSE_BOOST
