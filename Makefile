SRC_COMM=$(wildcard src/*.cpp) $(wildcard src/io/*.cpp)
OBJS_COMM=$(SRC_COMM:%.cpp=%.o)

SRC=$(wildcard *.cpp)
EXE=HBT HBTdouble HBTboost

default: HBTboost
include Makefile.inc

# $(EXE): $(OBJS_COMM)

HBTdouble: CXXFLAGS+=-DHBT_REAL8 -DHBT_INT8 
HBTboost:CXX=mpic++
HBTboost:src/config_parser.o src/io/snapshot_io.o
HBTboost: LDLIBS+=-lboost_serialization -lboost_mpi
HBTdouble: HBT.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

depend:
	makedepend --$(CXXFLAGS)-- -Y $(SRC)
	
synccosma: clean
	rsync -avzL $(shell pwd)/ jvbq85@cosma-a:data/HBT2/code
	
synccosmalocal: clean
	rsync -avzL -e "ssh -p 4800" $(shell pwd)/ jvbq85@localhost:data/HBT2/code	
