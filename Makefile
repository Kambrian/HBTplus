SRC=$(wildcard *.cpp)
EXE=HBT HBTdouble
OBJS_COMM=mymath.o config_parser.o simulation_io/snapshot.o simulation_io/halo.o simulation_io/subhalo.o gravity/tree.o

default: HBT
include Makefile.inc

HBTdouble: CXXFLAGS+=-DHBT_REAL8 -DHBT_INT8 

HBT HBTdouble: HBT.o $(OBJS_COMM)
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

depend:
	makedepend --$(CXXFLAGS)-- -Y $(SRC)
	
synccosma: clean
	rsync -avzL $(shell pwd)/ jvbq85@cosma-a:data/HBT2/code
	
synccosmalocal: clean
	rsync -avzL -e "ssh -p 4800" $(shell pwd)/ jvbq85@localhost:data/HBT2/code	
