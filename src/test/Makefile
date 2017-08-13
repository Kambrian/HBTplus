SRC=$(wildcard *.cpp)
EXE=$(basename $(SRC))
default: $(EXE)
test_parse:  ../config_parser.o

OUTPUTS=$(wildcard *.hdf5)

HBTDIR=../..
include $(HBTDIR)/Makefile.inc
test_omptask:CXXFLAGS+=-O0

