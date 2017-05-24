SRC_Jing=$(wildcard src/io/jing/*.cpp) $(wildcard src/io/jing/*.f90)
SRC_COMM=$(wildcard src/*.cpp) $(wildcard src/io/*.cpp) $(SRC_Jing)
OBJS_COMM=$(patsubst %.f90,%.f.o, $(SRC_COMM:%.cpp=%.o))

SRC=$(wildcard *.cpp)
EXE_HBT=HBT HBTdouble  HBT_majormerger_test  HBTi8 HBT.apostle HBT.apostle_thermal
EXE=$(EXE_HBT) FoF

default: HBT
include Makefile.inc

echo:
	@echo $(OBJS_COMM)

$(EXE): $(OBJS_COMM)

HBT.apostle HBT.apostle_thermal: CXXFLAGS+=-DHBT_INT8 -DUNSIGNED_LONG_ID_OUTPUT
HBT.apostle_thermal: CXXFLAGS+=-DUNBIND_WITH_THERMAL_ENERGY

HBTi8: CXXFLAGS+=-DHBT_INT8
HBTdouble: CXXFLAGS+=-DHBT_REAL8 -DHBT_INT8 
HBT_majormerger_test: CXXFLAGS+=-DMAJOR_MERGER_PATCH #-DALLOW_BINARY_SYSTEM
$(EXE_HBT): HBT.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

FoF: OMPFLAG=
FoF: CXXFLAGS+=-DDM_ONLY
	
depend:
	makedepend --$(CXXFLAGS)-- -Y $(SRC) $(SRC_COMM)

#custom command, not needed by a general user
-include .Makefile_sync_hydro.inc

# DO NOT DELETE

debug.o: src/datatypes.h src/config_parser.h src/datatypes.h src/snapshot.h
debug.o: src/mymath.h src/config_parser.h src/snapshot_number.h src/hash.h
debug.o: src/hash.tpp src/halo.h src/snapshot.h src/subhalo.h src/halo.h
debug.o: src/mymath.h
HBT.o: src/datatypes.h src/config_parser.h src/datatypes.h src/snapshot.h
HBT.o: src/mymath.h src/config_parser.h src/snapshot_number.h src/hash.h
HBT.o: src/hash.tpp src/halo.h src/snapshot.h src/subhalo.h src/halo.h
HBT.o: src/mymath.h
src/config_parser.o: src/config_parser.h src/datatypes.h
src/gravity_tree.o: src/mymath.h src/datatypes.h src/config_parser.h
src/gravity_tree.o: src/gravity_tree.h src/snapshot.h
src/halo.o: src/mymath.h src/datatypes.h src/config_parser.h src/halo.h
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
src/io/snapshot_io.o: src/snapshot.h src/datatypes.h src/mymath.h
src/io/snapshot_io.o: src/config_parser.h src/snapshot_number.h src/hash.h
src/io/snapshot_io.o: src/hash.tpp src/mymath.h
src/io/subhalo_io.o: src/datatypes.h src/snapshot_number.h src/datatypes.h
src/io/subhalo_io.o: src/config_parser.h src/subhalo.h src/snapshot_number.h
src/io/subhalo_io.o: src/halo.h
