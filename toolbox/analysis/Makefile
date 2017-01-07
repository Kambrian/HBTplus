SRC=$(wildcard *.cpp)
EXE=$(basename $(SRC))
#SRC_Jing=$(wildcard ../../src/io/jing/*.cpp) $(wildcard ../../src/io/jing/*.f90)
SRC_COMM=$(wildcard ../../src/*.cpp) $(wildcard ../../src/io/*.cpp) #$(SRC_Jing)
OBJS_COMM=$(patsubst %.f90,%.f.o, $(SRC_COMM:%.cpp=%.o))

default: $(EXE)

$(EXE): $(OBJS_COMM)

$(EXE): CXXFLAGS+=-I.. -DHBT_INT8 #-DHBT_REAL8

include ../../Makefile.inc
-include .Makefile_sync.inc
