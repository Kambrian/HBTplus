SRC_COMM=$(wildcard src/*.cpp) $(wildcard src/io/*.cpp)
OBJS_COMM=$(SRC_COMM:%.cpp=%.o)

SRC=$(wildcard *.cpp)
EXE_HBT=HBT HBTdouble  HBT_majormerger_test HBTi8
EXE=$(EXE_HBT) 
#EXE+=debug

default: HBT
include Makefile.inc

$(EXE): $(OBJS_COMM)

# debug: CXXFLAGS+=-DHBT_INT8 -DHBT_REAL8
HBTdouble: CXXFLAGS+=-DHBT_INT8 -DHBT_REAL8 
HBTi8: CXXFLAGS+=-DHBT_INT8
HBT_majormerger_test: CXXFLAGS+=-DMAJOR_MERGER_PATCH #-DALLOW_BINARY_SYSTEM
$(EXE_HBT): HBT.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

depend:
	makedepend --$(CXXFLAGS)-- -Y $(SRC) $(SRC_COMM)
	
synccosma: clean
	rsync -avzL $(shell pwd)/ jvbq85@cosma-c:data/HBT2/code
	
synccosmalocal: clean
	rsync -avzL -e "ssh -p 4800" $(shell pwd)/ jvbq85@localhost:data/HBT2/code	

# DO NOT DELETE

debug.o: src/mpi_wrapper.h src/datatypes.h src/datatypes.h
debug.o: src/config_parser.h src/mpi_wrapper.h src/snapshot.h src/mymath.h
debug.o: src/config_parser.h src/snapshot_number.h src/hash.h src/hash.tpp
debug.o: src/halo.h src/snapshot.h src/subhalo.h src/halo.h src/mymath.h
HBT.o: src/mpi_wrapper.h src/datatypes.h src/datatypes.h src/config_parser.h
HBT.o: src/mpi_wrapper.h src/snapshot.h src/mymath.h src/config_parser.h
HBT.o: src/snapshot_number.h src/hash.h src/hash.tpp src/halo.h
HBT.o: src/snapshot.h src/subhalo.h src/halo.h src/mymath.h
src/config_parser.o: src/config_parser.h
src/gravity_tree.o: src/mymath.h src/datatypes.h src/config_parser.h
src/gravity_tree.o: src/gravity_tree.h src/snapshot.h
src/halo.o: src/mpi_wrapper.h src/mymath.h src/datatypes.h
src/halo.o: src/config_parser.h src/halo.h
src/mpi_wrapper.o: src/mpi_wrapper.h
src/mymath.o: src/mymath.h src/datatypes.h src/config_parser.h
src/snapshot.o: src/snapshot.h src/mymath.h src/datatypes.h
src/snapshot.o: src/config_parser.h
src/subhalo.o: src/datatypes.h src/snapshot_number.h src/config_parser.h
src/subhalo.o: src/subhalo.h src/halo.h src/gravity_tree.h src/snapshot.h
src/subhalo_tracking.o: src/datatypes.h src/snapshot_number.h
src/subhalo_tracking.o: src/config_parser.h src/subhalo.h src/halo.h
src/subhalo_unbind.o: src/datatypes.h src/snapshot_number.h
src/subhalo_unbind.o: src/config_parser.h src/subhalo.h src/halo.h
src/subhalo_unbind.o: src/gravity_tree.h src/snapshot.h
src/io/halo_io.o: src/mymath.h src/datatypes.h src/config_parser.h src/halo.h
src/io/halo_io.o: src/snapshot_number.h src/snapshot.h
src/io/snapshot_io.o: src/mpi_wrapper.h src/datatypes.h src/snapshot.h
src/io/snapshot_io.o: src/mymath.h src/config_parser.h src/snapshot_number.h
src/io/snapshot_io.o: src/hash.h src/hash.tpp src/mpi_wrapper.h src/mymath.h
src/io/subhalo_io.o: src/mpi_wrapper.h src/datatypes.h src/datatypes.h
src/io/subhalo_io.o: src/snapshot_number.h src/config_parser.h src/subhalo.h
src/io/subhalo_io.o: src/snapshot_number.h src/halo.h
