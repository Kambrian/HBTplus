SRC=$(wildcard *.cpp)
EXE=$(basename $(SRC))
OBJS_COMM=../mymath.o ../config_parser.o ../simulation_io/snapshot.o ../simulation_io/halo.o ../simulation_io/subhalo.o ../gravity/tree.o

default: $(EXE)

$(EXE): $(OBJS_COMM)

include ../Makefile.inc

