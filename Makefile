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
	makedepend --$(CXXFLAGS)-- -Y $(SRC)
	
synccosma: clean
	rsync -avzL $(shell pwd)/ jvbq85@cosma-c:data/HBT2/code
	
synccosmalocal: clean
	rsync -avzL -e "ssh -p 4800" $(shell pwd)/ jvbq85@localhost:data/HBT2/code	
