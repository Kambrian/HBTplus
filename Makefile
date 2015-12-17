SRC_COMM=$(wildcard src/*.cpp) $(wildcard src/io/*.cpp)
OBJS_COMM=$(SRC_COMM:%.cpp=%.o)

SRC=$(wildcard *.cpp)
EXE_VARIANT=HBTdouble  HBT_majormerger_test HBTmajormerger
EXE=HBT $(EXE_VARIANT)

default: HBT
include Makefile.inc

$(EXE): $(OBJS_COMM)

HBTdouble: CXXFLAGS+=-DHBT_REAL8 -DHBT_INT8 
HBTmajormerger: CXXFLAGS+=-DALLOW_BINARY_SYSTEM
HBT_majormerger_test: CXXFLAGS+=-DMAJOR_MERGER_PATCH -DALLOW_BINARY_SYSTEM
$(EXE_VARIANT): HBT.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@

depend:
	makedepend --$(CXXFLAGS)-- -Y $(SRC)
	
synccosma: clean
	rsync -avzL $(shell pwd)/ jvbq85@cosma-a:data/HBT2/code
	
synccosmalocal: clean
	rsync -avzL -e "ssh -p 4800" $(shell pwd)/ jvbq85@localhost:data/HBT2/code	
