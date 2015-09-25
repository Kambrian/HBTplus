include Makefile.inc
SRC=$(wildcard *.cpp)
EXE=HBT

HBT: mymath.o config_parser.o simulation_io/snapshot.o simulation_io/halo.o simulation_io/subhalo.o

depend:
	makedepend --$(CXXFLAGS)-- -Y $(SRC)
	
synccosma: clean
	rsync -avzL $(shell pwd)/ jvbq85@cosma-a:data/HBT2/code# DO NOT DELETE
