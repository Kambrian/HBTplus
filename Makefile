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

./src/config_parser.o: ./src/config_parser.h ./src/datatypes.h ./src/mymath.h
./src/config_parser.o: ./src/mysort.h
./src/subhalo_tracking.o: ./src/datatypes.h ./src/snapshot_number.h
./src/subhalo_tracking.o: ./src/config_parser.h ./src/subhalo.h ./src/halo.h
./src/subhalo_tracking.o: ./src/snapshot.h ./src/mymath.h ./src/mysort.h
./src/subhalo_tracking.o: ./src/hash.h ./src/hash.tpp ./src/hdf_wrapper.h
./src/snapshot.o: ./src/snapshot.h ./src/datatypes.h ./src/mymath.h
./src/snapshot.o: ./src/mysort.h ./src/config_parser.h
./src/snapshot.o: ./src/snapshot_number.h ./src/hash.h ./src/hash.tpp
./src/linkedlist.o: ./src/linkedlist.h ./src/mymath.h ./src/mysort.h
./src/linkedlist.o: ./src/datatypes.h ./src/linkedlist_base.h
./src/linkedlist.o: ./src/snapshot.h ./src/config_parser.h
./src/linkedlist.o: ./src/snapshot_number.h ./src/hash.h ./src/hash.tpp
./src/hdf_wrapper.o: ./src/hdf_wrapper.h
./src/halo.o: ./src/mymath.h ./src/mysort.h ./src/datatypes.h ./src/halo.h
./src/halo.o: ./src/snapshot_number.h ./src/config_parser.h ./src/snapshot.h
./src/halo.o: ./src/hash.h ./src/hash.tpp
./src/subhalo_merge.o: ./src/datatypes.h ./src/snapshot_number.h
./src/subhalo_merge.o: ./src/config_parser.h ./src/subhalo.h ./src/halo.h
./src/subhalo_merge.o: ./src/snapshot.h ./src/mymath.h ./src/mysort.h
./src/subhalo_merge.o: ./src/hash.h ./src/hash.tpp ./src/hdf_wrapper.h
./src/geometric_tree.o: ./src/mymath.h ./src/mysort.h ./src/datatypes.h
./src/geometric_tree.o: ./src/config_parser.h ./src/geometric_tree.h
./src/geometric_tree.o: ./src/oct_tree.h ./src/snapshot.h
./src/geometric_tree.o: ./src/snapshot_number.h ./src/hash.h ./src/hash.tpp
./src/geometric_tree.o: ./src/oct_tree.tpp
./src/subhalo.o: ./src/datatypes.h ./src/snapshot_number.h
./src/subhalo.o: ./src/config_parser.h ./src/subhalo.h ./src/halo.h
./src/subhalo.o: ./src/snapshot.h ./src/mymath.h ./src/mysort.h ./src/hash.h
./src/subhalo.o: ./src/hash.tpp ./src/hdf_wrapper.h ./src/gravity_tree.h
./src/subhalo.o: ./src/oct_tree.h ./src/oct_tree.tpp
./src/gravity_tree.o: ./src/mymath.h ./src/mysort.h ./src/datatypes.h
./src/gravity_tree.o: ./src/config_parser.h ./src/gravity_tree.h
./src/gravity_tree.o: ./src/oct_tree.h ./src/snapshot.h
./src/gravity_tree.o: ./src/snapshot_number.h ./src/hash.h ./src/hash.tpp
./src/gravity_tree.o: ./src/oct_tree.tpp
./src/mymath.o: ./src/mymath.h ./src/mysort.h ./src/datatypes.h
./src/subhalo_unbind.o: ./src/datatypes.h ./src/snapshot_number.h
./src/subhalo_unbind.o: ./src/config_parser.h ./src/subhalo.h ./src/halo.h
./src/subhalo_unbind.o: ./src/snapshot.h ./src/mymath.h ./src/mysort.h
./src/subhalo_unbind.o: ./src/hash.h ./src/hash.tpp ./src/hdf_wrapper.h
./src/subhalo_unbind.o: ./src/gravity_tree.h ./src/oct_tree.h
./src/subhalo_unbind.o: ./src/oct_tree.tpp
./src/linkedlist_base.o: ./src/mymath.h ./src/mysort.h ./src/datatypes.h
./src/linkedlist_base.o: ./src/linkedlist_base.h ./src/snapshot.h
./src/linkedlist_base.o: ./src/config_parser.h ./src/snapshot_number.h
./src/linkedlist_base.o: ./src/hash.h ./src/hash.tpp
./src/fof_builder.o: ./src/fof_builder.h ./src/geometric_tree.h
./src/fof_builder.o: ./src/oct_tree.h ./src/datatypes.h ./src/snapshot.h
./src/fof_builder.o: ./src/mymath.h ./src/mysort.h ./src/config_parser.h
./src/fof_builder.o: ./src/snapshot_number.h ./src/hash.h ./src/hash.tpp
./src/fof_builder.o: ./src/oct_tree.tpp
./src/io/custom_io.o: ./src/mymath.h ./src/mysort.h ./src/datatypes.h
./src/io/custom_io.o: ./src/io/custom_io.h ./src/datatypes.h ./src/halo.h
./src/io/custom_io.o: ./src/snapshot_number.h ./src/config_parser.h
./src/io/custom_io.o: ./src/snapshot.h ./src/mymath.h ./src/hash.h
./src/io/custom_io.o: ./src/hash.tpp
./src/io/hbt_group_io.o: ./src/mymath.h ./src/mysort.h ./src/datatypes.h
./src/io/hbt_group_io.o: ./src/io/hbt_group_io.h ./src/datatypes.h
./src/io/hbt_group_io.o: ./src/halo.h ./src/snapshot_number.h
./src/io/hbt_group_io.o: ./src/config_parser.h ./src/snapshot.h
./src/io/hbt_group_io.o: ./src/mymath.h ./src/hash.h ./src/hash.tpp
./src/io/hbt_group_io.o: ./src/hdf_wrapper.h
./src/io/snapshot_io.o: ./src/snapshot.h ./src/datatypes.h ./src/mymath.h
./src/io/snapshot_io.o: ./src/mysort.h ./src/config_parser.h
./src/io/snapshot_io.o: ./src/snapshot_number.h ./src/hash.h ./src/hash.tpp
./src/io/snapshot_io.o: ./src/mymath.h ./src/io/gadget_io.h
./src/io/snapshot_io.o: ./src/io/apostle_io.h ./src/hdf_wrapper.h
./src/io/snapshot_io.o: ./src/halo.h ./src/snapshot.h
./src/io/gadget_io.o: ./src/snapshot.h ./src/datatypes.h ./src/mymath.h
./src/io/gadget_io.o: ./src/mysort.h ./src/config_parser.h
./src/io/gadget_io.o: ./src/snapshot_number.h ./src/hash.h ./src/hash.tpp
./src/io/gadget_io.o: ./src/mymath.h ./src/io/gadget_io.h
./src/io/apostle_io.o: ./src/snapshot.h ./src/datatypes.h ./src/mymath.h
./src/io/apostle_io.o: ./src/mysort.h ./src/config_parser.h
./src/io/apostle_io.o: ./src/snapshot_number.h ./src/hash.h ./src/hash.tpp
./src/io/apostle_io.o: ./src/mymath.h ./src/hdf_wrapper.h
./src/io/apostle_io.o: ./src/io/apostle_io.h ./src/halo.h ./src/snapshot.h
./src/io/subhalo_io.o: ./src/datatypes.h ./src/snapshot_number.h
./src/io/subhalo_io.o: ./src/datatypes.h ./src/config_parser.h
./src/io/subhalo_io.o: ./src/subhalo.h ./src/snapshot_number.h ./src/halo.h
./src/io/subhalo_io.o: ./src/snapshot.h ./src/mymath.h ./src/mysort.h
./src/io/subhalo_io.o: ./src/hash.h ./src/hash.tpp ./src/hdf_wrapper.h
./src/io/subhalo_io.o: ./src/hdf_wrapper.h
./src/io/gadget_group_io.o: ./src/mymath.h ./src/mysort.h ./src/datatypes.h
./src/io/gadget_group_io.o: ./src/io/gadget_group_io.h ./src/halo.h
./src/io/gadget_group_io.o: ./src/snapshot_number.h ./src/config_parser.h
./src/io/gadget_group_io.o: ./src/snapshot.h ./src/mymath.h ./src/hash.h
./src/io/gadget_group_io.o: ./src/hash.tpp
./src/io/halo_io.o: ./src/mymath.h ./src/mysort.h ./src/datatypes.h
./src/io/halo_io.o: ./src/halo.h ./src/snapshot_number.h
./src/io/halo_io.o: ./src/config_parser.h ./src/snapshot.h ./src/mymath.h
./src/io/halo_io.o: ./src/hash.h ./src/hash.tpp ./src/io/gadget_group_io.h
./src/io/halo_io.o: ./src/io/apostle_io.h ./src/hdf_wrapper.h
./src/io/halo_io.o: ./src/io/hbt_group_io.h ./src/datatypes.h
./src/io/halo_io.o: ./src/io/custom_io.h
./src/io/subhalo_io_multiple.o: ./src/datatypes.h ./src/snapshot_number.h
./src/io/subhalo_io_multiple.o: ./src/datatypes.h ./src/config_parser.h
./src/io/subhalo_io_multiple.o: ./src/subhalo.h ./src/snapshot_number.h
./src/io/subhalo_io_multiple.o: ./src/halo.h ./src/snapshot.h ./src/mymath.h
./src/io/subhalo_io_multiple.o: ./src/mysort.h ./src/hash.h ./src/hash.tpp
./src/io/subhalo_io_multiple.o: ./src/hdf_wrapper.h
