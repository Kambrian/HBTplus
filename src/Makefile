SRC=$(wildcard *.cpp)
OBJS=$(SRC:%.cpp=%.o)
TEST=$(SRC:%.cpp=%.test)

gravity_tree.test: config_parser.o mymath.o halo.o snapshot.o io/halo_io.o io/snapshot_io.o

$(TEST): HDFLIB=

include ../Makefile.inc
	
