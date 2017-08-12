EXE_HBT=HBT HBTdouble  HBT_majormerger_test  HBTi8 HBT.apostle HBT.apostle_thermal HBT.nostrip
EXE_FOF=FoF FoF.ll FoFdebug FoFdebug2
EXE=$(EXE_HBT) $(EXE_FOF)

default: HBT

HBTDIR=.
include $(HBTDIR)/Makefile.inc

$(EXE): $(OBJS_COMM)

HBT.apostle HBT.apostle_thermal: CXXFLAGS+=-DHBT_INT8 -DUNSIGNED_LONG_ID_OUTPUT
HBT.apostle_thermal: CXXFLAGS+=-DUNBIND_WITH_THERMAL_ENERGY

HBTi8: CXXFLAGS+=-DHBT_INT8
HBTdouble: CXXFLAGS+=-DHBT_REAL8 -DHBT_INT8 
HBT_majormerger_test: CXXFLAGS+=-DMAJOR_MERGER_PATCH #-DALLOW_BINARY_SYSTEM
HBT.nostrip: CXXFLAGS+=-DNO_STRIPPING -DHBT_INT8 #keep tracking and do not remove unbound particles
$(EXE_HBT): HBT.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

FoF.ll: CXXFLAGS+=-DFOF_METHOD=2
$(EXE_FOF): OMPFLAG=
$(EXE_FOF): CXXFLAGS+=-DDM_ONLY -DHBT_INT8
# $(EXE_FOF): FoF.o
# 	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@
# 	rm FoF.o


depend:
	makedepend --$(CXXFLAGS)-- -Y $(SRC) $(SRC_COMM)

#custom command, not needed by a general user
-include .Makefile_sync_hydro.inc

# DO NOT DELETE

FoF.o: src/datatypes.h src/config_parser.h src/datatypes.h src/snapshot.h
FoF.o: src/mymath.h src/config_parser.h src/snapshot_number.h src/hash.h
FoF.o: src/hash.tpp src/io/hbt_group_io.h src/datatypes.h src/halo.h
FoF.o: src/snapshot.h src/linkedlist.h src/fof_builder.h src/geometric_tree.h
FoF.o: src/oct_tree.h src/oct_tree.tpp
HBT.o: src/datatypes.h src/config_parser.h src/datatypes.h src/snapshot.h
HBT.o: src/mymath.h src/config_parser.h src/snapshot_number.h src/hash.h
HBT.o: src/hash.tpp src/halo.h src/snapshot.h src/subhalo.h src/halo.h
HBT.o: src/hdf_wrapper.h src/mymath.h
FoFdebug2.o: src/datatypes.h src/config_parser.h src/datatypes.h
FoFdebug2.o: src/snapshot.h src/mymath.h src/config_parser.h
FoFdebug2.o: src/snapshot_number.h src/hash.h src/hash.tpp src/halo.h
FoFdebug2.o: src/snapshot.h src/subhalo.h src/halo.h src/hdf_wrapper.h
FoFdebug2.o: src/mymath.h src/io/hbt_group_io.h src/datatypes.h src/halo.h
FoFdebug2.o: src/gravity_tree.h src/oct_tree.h src/oct_tree.tpp
FoFdebug2.o: src/linkedlist.h src/fof_builder.h src/geometric_tree.h
FoFdebug.o: src/datatypes.h src/config_parser.h src/datatypes.h
FoFdebug.o: src/snapshot.h src/mymath.h src/config_parser.h
FoFdebug.o: src/snapshot_number.h src/hash.h src/hash.tpp src/halo.h
FoFdebug.o: src/snapshot.h src/subhalo.h src/halo.h src/hdf_wrapper.h
FoFdebug.o: src/mymath.h src/io/hbt_group_io.h src/datatypes.h src/halo.h
FoFdebug.o: src/gravity_tree.h src/oct_tree.h src/oct_tree.tpp
FoFdebug.o: src/linkedlist.h src/fof_builder.h src/geometric_tree.h
HBTdebug.o: src/datatypes.h src/config_parser.h src/datatypes.h
HBTdebug.o: src/snapshot.h src/mymath.h src/config_parser.h
HBTdebug.o: src/snapshot_number.h src/hash.h src/hash.tpp src/halo.h
HBTdebug.o: src/snapshot.h src/subhalo.h src/halo.h src/hdf_wrapper.h
HBTdebug.o: src/mymath.h
src/config_parser.o: src/config_parser.h src/mymath.h src/datatypes.h
src/subhalo_tracking.o: src/datatypes.h src/snapshot_number.h
src/subhalo_tracking.o: src/config_parser.h src/subhalo.h src/halo.h
src/subhalo_tracking.o: src/hdf_wrapper.h
src/snapshot.o: src/snapshot.h src/mymath.h src/datatypes.h
src/linkedlist.o: src/mymath.h src/datatypes.h src/linkedlist.h
src/linkedlist.o: src/snapshot.h
src/hdf_wrapper.o: src/hdf_wrapper.h
src/subhalo_merge.o: src/datatypes.h src/snapshot_number.h
src/subhalo_merge.o: src/config_parser.h src/subhalo.h src/halo.h
src/subhalo_merge.o: src/hdf_wrapper.h
src/geometric_tree.o: src/mymath.h src/datatypes.h src/config_parser.h
src/geometric_tree.o: src/geometric_tree.h src/oct_tree.h src/snapshot.h
src/geometric_tree.o: src/oct_tree.tpp
src/gravity_tree.o: src/mymath.h src/datatypes.h src/config_parser.h
src/gravity_tree.o: src/gravity_tree.h src/oct_tree.h src/snapshot.h
src/gravity_tree.o: src/oct_tree.tpp
src/subhalo.o: src/datatypes.h src/snapshot_number.h src/config_parser.h
src/subhalo.o: src/subhalo.h src/halo.h src/hdf_wrapper.h src/gravity_tree.h
src/subhalo.o: src/oct_tree.h src/snapshot.h src/oct_tree.tpp src/mymath.h
src/mymath.o: src/mymath.h src/datatypes.h
src/subhalo_unbind.o: src/datatypes.h src/snapshot_number.h
src/subhalo_unbind.o: src/config_parser.h src/subhalo.h src/halo.h
src/subhalo_unbind.o: src/hdf_wrapper.h src/gravity_tree.h src/oct_tree.h
src/subhalo_unbind.o: src/snapshot.h src/oct_tree.tpp src/mymath.h
src/halo.o: src/mymath.h src/datatypes.h src/halo.h
src/linkedlist_parallel.o: src/linkedlist_parallel.h src/mymath.h
src/linkedlist_parallel.o: src/datatypes.h src/linkedlist.h src/snapshot.h
src/fof_builder.o: src/fof_builder.h src/geometric_tree.h src/oct_tree.h
src/fof_builder.o: src/datatypes.h src/snapshot.h src/oct_tree.tpp
src/fof_builder.o: src/mymath.h src/config_parser.h
src/io/custom_io.o: src/mymath.h src/datatypes.h src/io/custom_io.h
src/io/custom_io.o: src/datatypes.h src/halo.h src/snapshot_number.h
src/io/custom_io.o: src/config_parser.h src/snapshot.h
src/io/hbt_group_io.o: src/mymath.h src/datatypes.h src/io/hbt_group_io.h
src/io/hbt_group_io.o: src/datatypes.h src/halo.h src/snapshot_number.h
src/io/hbt_group_io.o: src/config_parser.h src/snapshot.h src/hdf_wrapper.h
src/io/snapshot_io.o: src/snapshot.h src/datatypes.h src/mymath.h
src/io/snapshot_io.o: src/config_parser.h src/snapshot_number.h src/hash.h
src/io/snapshot_io.o: src/hash.tpp src/mymath.h src/io/gadget_io.h
src/io/snapshot_io.o: src/io/apostle_io.h src/hdf_wrapper.h src/halo.h
src/io/snapshot_io.o: src/snapshot.h src/io/jing/jing_io.h src/halo.h
src/io/gadget_io.o: src/snapshot.h src/datatypes.h src/mymath.h
src/io/gadget_io.o: src/config_parser.h src/snapshot_number.h src/hash.h
src/io/gadget_io.o: src/hash.tpp src/mymath.h src/io/gadget_io.h
src/io/apostle_io.o: src/snapshot.h src/datatypes.h src/mymath.h
src/io/apostle_io.o: src/config_parser.h src/snapshot_number.h src/hash.h
src/io/apostle_io.o: src/hash.tpp src/mymath.h src/hdf_wrapper.h
src/io/apostle_io.o: src/io/apostle_io.h src/halo.h src/snapshot.h
src/io/subhalo_io.o: src/datatypes.h src/snapshot_number.h src/datatypes.h
src/io/subhalo_io.o: src/config_parser.h src/subhalo.h src/snapshot_number.h
src/io/subhalo_io.o: src/halo.h src/hdf_wrapper.h src/hdf_wrapper.h
src/io/gadget_group_io.o: src/mymath.h src/datatypes.h
src/io/gadget_group_io.o: src/io/gadget_group_io.h src/halo.h
src/io/gadget_group_io.o: src/snapshot_number.h src/config_parser.h
src/io/gadget_group_io.o: src/snapshot.h
src/io/halo_io.o: src/mymath.h src/datatypes.h src/halo.h
src/io/halo_io.o: src/snapshot_number.h src/config_parser.h src/snapshot.h
src/io/halo_io.o: src/io/gadget_group_io.h src/io/apostle_io.h
src/io/halo_io.o: src/hdf_wrapper.h src/io/jing/jing_io.h src/halo.h
src/io/halo_io.o: src/io/hbt_group_io.h src/datatypes.h src/io/custom_io.h
src/io/subhalo_io_multiple.o: src/datatypes.h src/snapshot_number.h
src/io/subhalo_io_multiple.o: src/datatypes.h src/config_parser.h
src/io/subhalo_io_multiple.o: src/subhalo.h src/snapshot_number.h src/halo.h
src/io/subhalo_io_multiple.o: src/hdf_wrapper.h
src/io/jing/jing_io.o: src/snapshot.h src/datatypes.h src/mymath.h
src/io/jing/jing_io.o: src/config_parser.h src/snapshot_number.h src/hash.h
src/io/jing/jing_io.o: src/hash.tpp src/io/jing/jing_io.h src/halo.h
src/io/jing/jing_io.o: src/io/jing/fortread.h
