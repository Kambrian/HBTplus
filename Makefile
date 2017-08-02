SRC_COMM=$(wildcard src/*.cpp) $(wildcard src/io/*.cpp) 
OBJS_COMM=$(SRC_COMM:%.cpp=%.o)

SRC=$(wildcard *.cpp)
EXE_HBT=HBT HBTdouble  HBT_majormerger_test  HBTi8 HBT.apostle HBT.apostle_thermal HBT.nostrip
EXE=$(EXE_HBT)
# EXE+=debug

default: HBT
include Makefile.inc

$(EXE): $(OBJS_COMM)

HBT.apostle HBT.apostle_thermal: CXXFLAGS+=-DHBT_INT8 -DUNSIGNED_LONG_ID_OUTPUT
HBT.apostle_thermal: CXXFLAGS+=-DUNBIND_WITH_THERMAL_ENERGY
# debug: CXXFLAGS+=-DHBT_INT8 -DHBT_REAL8

HBTdouble: CXXFLAGS+=-DHBT_INT8 -DHBT_REAL8 
HBTi8: CXXFLAGS+=-DHBT_INT8
HBT_majormerger_test: CXXFLAGS+=-DMAJOR_MERGER_PATCH #-DALLOW_BINARY_SYSTEM
HBT.nostrip: CXXFLAGS+=-DNO_STRIPPING -DHBT_INT8 #track without unbinding.
$(EXE_HBT): HBT.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

depend:
	makedepend --$(CXXFLAGS)-- -Y $(SRC) $(SRC_COMM)
	
#custom command, not needed by a general user
-include .Makefile_sync_mpi.inc
# DO NOT DELETE

HBT.o: src/mpi_wrapper.h src/datatypes.h src/mymath.h src/datatypes.h
HBT.o: src/config_parser.h src/mpi_wrapper.h src/snapshot.h
HBT.o: src/config_parser.h src/snapshot_number.h src/hash.h src/hash.tpp
HBT.o: src/halo.h src/snapshot.h src/subhalo.h src/halo.h src/hdf_wrapper.h
HBT.o: src/mymath.h src/particle_exchanger.h src/halo_particle_iterator.h
HBT.o: src/subhalo.h src/hash_remote.tpp
FoFdebug2.o: src/datatypes.h src/config_parser.h src/datatypes.h
FoFdebug2.o: src/mpi_wrapper.h src/snapshot.h src/mymath.h
FoFdebug2.o: src/config_parser.h src/snapshot_number.h src/hash.h
FoFdebug2.o: src/hash.tpp src/halo.h src/snapshot.h src/subhalo.h src/halo.h
FoFdebug2.o: src/hdf_wrapper.h src/mymath.h src/gravity_tree.h src/oct_tree.h
FoFdebug2.o: src/oct_tree.tpp src/linkedlist.h
FoFdebug.o: src/datatypes.h src/config_parser.h src/datatypes.h
FoFdebug.o: src/mpi_wrapper.h src/snapshot.h src/mymath.h src/config_parser.h
FoFdebug.o: src/snapshot_number.h src/hash.h src/hash.tpp src/halo.h
FoFdebug.o: src/snapshot.h src/subhalo.h src/halo.h src/hdf_wrapper.h
FoFdebug.o: src/mymath.h src/gravity_tree.h src/oct_tree.h src/oct_tree.tpp
FoFdebug.o: src/linkedlist.h
HBTdebug.o: src/mpi_wrapper.h src/datatypes.h src/mymath.h src/datatypes.h
HBTdebug.o: src/config_parser.h src/mpi_wrapper.h src/snapshot.h
HBTdebug.o: src/config_parser.h src/snapshot_number.h src/hash.h src/hash.tpp
HBTdebug.o: src/halo.h src/snapshot.h src/subhalo.h src/halo.h
HBTdebug.o: src/hdf_wrapper.h src/mymath.h src/particle_exchanger.h
HBTdebug.o: src/halo_particle_iterator.h src/subhalo.h src/hash_remote.tpp
src/subhalo_unbind.o: src/datatypes.h src/snapshot_number.h
src/subhalo_unbind.o: src/config_parser.h src/subhalo.h src/gravity_tree.h
src/subhalo_unbind.o: src/oct_tree.h src/snapshot.h src/oct_tree.tpp
src/subhalo_unbind.o: src/mymath.h
src/particle_exchanger.o: src/snapshot.h src/particle_exchanger.h
src/particle_exchanger.o: src/datatypes.h src/mymath.h src/mpi_wrapper.h
src/particle_exchanger.o: src/halo_particle_iterator.h src/subhalo.h
src/particle_exchanger.o: src/hash_remote.tpp src/hash.h src/hash.tpp
src/snapshot.o: src/snapshot.h src/mymath.h src/datatypes.h
src/linkedlist.o: src/mymath.h src/datatypes.h src/linkedlist.h
src/linkedlist.o: src/snapshot.h
src/hdf_wrapper.o: src/hdf_wrapper.h
src/subhalo_merge.o: src/datatypes.h src/snapshot_number.h
src/subhalo_merge.o: src/config_parser.h src/subhalo.h
src/mpi_wrapper.o: src/mpi_wrapper.h
src/snapshot_exchanger.o: src/snapshot.h src/mymath.h src/datatypes.h
src/snapshot_exchanger.o: src/mpi_wrapper.h
src/geometric_tree.o: src/mymath.h src/datatypes.h src/config_parser.h
src/geometric_tree.o: src/geometric_tree.h src/oct_tree.h src/snapshot.h
src/geometric_tree.o: src/oct_tree.tpp
src/subhalo_tracking.o: src/datatypes.h src/snapshot_number.h
src/subhalo_tracking.o: src/config_parser.h src/subhalo.h
src/gravity_tree.o: src/mymath.h src/datatypes.h src/config_parser.h
src/gravity_tree.o: src/gravity_tree.h src/oct_tree.h src/snapshot.h
src/gravity_tree.o: src/oct_tree.tpp
src/subhalo.o: src/datatypes.h src/snapshot_number.h src/config_parser.h
src/subhalo.o: src/subhalo.h src/particle_exchanger.h src/mymath.h
src/subhalo.o: src/mpi_wrapper.h src/snapshot.h src/halo_particle_iterator.h
src/subhalo.o: src/hash_remote.tpp src/hash.h src/hash.tpp
src/mymath.o: src/mymath.h src/datatypes.h
src/halo.o: src/mpi_wrapper.h src/mymath.h src/datatypes.h src/halo.h
src/halo.o: src/particle_exchanger.h src/snapshot.h
src/halo.o: src/halo_particle_iterator.h src/subhalo.h src/hash_remote.tpp
src/halo.o: src/hash.h src/hash.tpp
src/linkedlist_parallel.o: src/linkedlist_parallel.h src/mymath.h
src/linkedlist_parallel.o: src/datatypes.h src/linkedlist.h src/snapshot.h
src/config_parser.o: src/config_parser.h
src/io/snapshot_io.o: src/mpi_wrapper.h src/datatypes.h src/mymath.h
src/io/snapshot_io.o: src/snapshot.h src/config_parser.h
src/io/snapshot_io.o: src/snapshot_number.h src/hash.h src/hash.tpp
src/io/snapshot_io.o: src/mpi_wrapper.h src/mymath.h src/io/gadget_io.h
src/io/snapshot_io.o: src/io/apostle_io.h src/hdf_wrapper.h src/halo.h
src/io/snapshot_io.o: src/snapshot.h
src/io/gadget_io.o: src/snapshot.h src/datatypes.h src/mymath.h
src/io/gadget_io.o: src/config_parser.h src/snapshot_number.h src/hash.h
src/io/gadget_io.o: src/hash.tpp src/mpi_wrapper.h src/mymath.h
src/io/gadget_io.o: src/io/gadget_io.h src/mpi_wrapper.h
src/io/apostle_io.o: src/snapshot.h src/datatypes.h src/mymath.h
src/io/apostle_io.o: src/config_parser.h src/snapshot_number.h src/hash.h
src/io/apostle_io.o: src/hash.tpp src/mpi_wrapper.h src/mymath.h
src/io/apostle_io.o: src/hdf_wrapper.h src/io/apostle_io.h src/halo.h
src/io/apostle_io.o: src/snapshot.h src/mpi_wrapper.h
src/io/apostle_io.o: src/halo_particle_iterator.h
src/io/subhalo_io.o: src/mpi_wrapper.h src/datatypes.h src/mymath.h
src/io/subhalo_io.o: src/datatypes.h src/snapshot_number.h
src/io/subhalo_io.o: src/config_parser.h src/subhalo.h src/snapshot_number.h
src/io/subhalo_io.o: src/halo.h src/hdf_wrapper.h
src/io/gadget_group_io.o: src/mymath.h src/halo.h src/datatypes.h
src/io/gadget_group_io.o: src/snapshot_number.h src/config_parser.h
src/io/gadget_group_io.o: src/snapshot.h src/io/gadget_group_io.h
src/io/gadget_group_io.o: src/mpi_wrapper.h src/mymath.h
src/io/halo_io.o: src/mymath.h src/halo.h src/datatypes.h
src/io/halo_io.o: src/snapshot_number.h src/config_parser.h src/snapshot.h
src/io/halo_io.o: src/io/gadget_group_io.h src/mpi_wrapper.h src/mymath.h
src/io/halo_io.o: src/io/apostle_io.h src/hdf_wrapper.h
