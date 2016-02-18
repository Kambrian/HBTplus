SRC=$(wildcard *.cpp)
EXE=$(basename $(SRC))
SRC_COMM=$(wildcard ../src/*.cpp) $(wildcard ../src/io/*.cpp)
OBJS_COMM=$(SRC_COMM:%.cpp=%.o)

default: $(EXE)

$(EXE): $(OBJS_COMM)

$(EXE): CXXFLAGS+=-I.. -DHBT_INT8 -DHBT_REAL8

include ../Makefile.inc

