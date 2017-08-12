CXX=g++

CXXFLAGS+=-std=c++11 $(HDFINC) $(GSLINC) $(OMPCFLAG) -O3 -g 
LDFLAGS+=-g
LDLIBS+=$(FTNLIB) $(HDFLIB) $(GSLLIB) $(OMPLIB) $(MATHLIB) 

LINK.o=$(CXX) $(LDFLAGS) $(TARGET_ARCH) #invoke c++ compiler for linking

#HDFINC=-I/usr/include/hdf5/serial/
#HDFLIB=-L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5_hl -lhdf5
HDFINC=
HDFLIB=-lhdf5_hl -lhdf5
#HDFINC+=-Wno-literal-suffix #to disable the annoying literal-suffix warning with the hdf header

###comment out the following lines if you do not have GSL (for eigenvalue decomposition of inertial tensors)
#GSLINC=
#GSLINC+=-DHAS_GSL 
#GSLLIB=-lgsl -lgslcblas
###comment up to here to disable GSL.

###for Dark Matter only simulations
CXXFLAGS+=-DDM_ONLY
###

##enable these if you want to obtain (thermal) binding energy for each particle
###save binding energy:
#CXXFLAGS+=-DSAVE_BINDING_ENERGY
###
###load thermal energy (also set SaveSubParticleProperties to 1 in config file if you want to save thermal energy)
#CXXFLAGS+=-DHAS_THERMAL_ENERGY
###

###for thermal unbinding:
#CXXFLAGS+=-DUNBIND_WITH_THERMAL_ENERGY
### alternatively, you can save both the non-thermal binding energy and the thermal energy and do it yourself. see above.

###To produce inclusive subhalos
#CXXFLAGS+=-DINCLUSIVE_MASS
###

ifeq ($(CXX), g++)
OMPCFLAG=-fopenmp
MATHLIB=-lm
endif
ifeq ($(CXX), icpc)
#new intel compilers accept only -qopenmp and -fopenmp; old intel only accept -openmp
OMPCFLAG=-qopenmp
MATHLIB=-limf
endif

OMPLIB=$(OMPCFLAG)

#CXXFLAGS+=-xHost #to get it runnable on cosma-e

include $(HBTDIR)/Makefile.jing_io

%.test.o: %.cpp
	$(CXX) $^ -DTEST_$(basename $<) $(CXXFLAGS) -c -o $@
%.test: %.test.o
	$(CXX) $^ $(LDFLAGS) $(LDLIBS) -o $@
$(TEST): CXXFLAGS+=-g -gdwarf-2 -O0

SUBDIRS:=$(wildcard */.)
SUBDIRSCLEAN=$(addsuffix clean,$(SUBDIRS))

clean: clean_curdir $(SUBDIRSCLEAN)

clean_curdir:
	rm -rf *.o *~ core a.out $(EXE) $(OUTPUTS) *.test

%clean: %
	$(MAKE) -C $< clean

depend_curdir:
	makedepend --$(CXXFLAGS)-- -Y $(SRC)
	
SRC_COMM=$(wildcard $(HBTDIR)/src/*.cpp) $(wildcard $(HBTDIR)/src/io/*.cpp) $(SRC_Jing)
OBJS_COMM=$(patsubst %.f90,%.f.o, $(SRC_COMM:%.cpp=%.o))

echo:
	@echo $(OBJS_COMM)