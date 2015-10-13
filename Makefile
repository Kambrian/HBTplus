SRC=$(wildcard *.cpp)
EXE=HBT
default: $(EXE)
include Makefile.inc

HBT: mymath.o config_parser.o simulation_io/snapshot.o simulation_io/halo.o simulation_io/subhalo.o gravity/tree.o

depend:
	makedepend --$(CXXFLAGS)-- -Y $(SRC)
	
synccosma: clean
	rsync -avzL $(shell pwd)/ jvbq85@cosma-a:data/HBT2/code
	
synccosmalocal: clean
	rsync -avzL -e "ssh -p 4800" $(shell pwd)/ jvbq85@localhost:data/HBT2/code	
