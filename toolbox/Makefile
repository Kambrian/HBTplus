SRC=$(wildcard *.cpp)
EXE=$(basename $(SRC))
default: $(EXE)

HBTDIR=..
include $(HBTDIR)/Makefile.inc

$(EXE): $(OBJS_COMM)
$(EXE): CXXFLAGS+=-I.. -DHBT_INT8 #-DHBT_REAL8




