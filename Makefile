SHELL      = /bin/sh

VPATH      = src/LAMMPS_STUB src/main src/mathutils src/stiffness_kernels src/surfaces src/unittests src/unittests/gtest-1.6.0/src/

CC         = g++
CCFLAGS    = -g -O0
#CCFLAGS    = -O0 -g `nc-config --cflags`
LINK       = g++
LINKFLAGS  = -g -O0 -pthread #`nc-config --libs` -lnetcdf_c++ 

INC        = -Isrc/LAMMPS_STUB -Isrc/main -Isrc/mathutils -Isrc/stiffness_kernels -Isrc/surfaces -Isrc/unittests/gtest-1.6.0 -Isrc/unittests/gtest-1.6.0/include
LIB        = 

SRC        = \
	memory.cpp \
	linearalgebra.cpp \
	gfmd_memory.cpp \
	gfmd_stiffness.cpp \
	table2d.cpp

TESTSRC = \
	memory.cpp \
	linearalgebra.cpp \
	boussinesq_cerruti.cpp \
	crystal_surface.cpp \
	dia100_surface.cpp \
	dia111_surface.cpp \
	fcc100_surface.cpp \
	fcc111_surface.cpp \
	sc100_surface.cpp \
	geometry.cpp \
	li_berger.cpp \
	pohrt_li.cpp \
	gtest-all.cc \
	test_crystal_surface.cpp \
	test_geometry.cpp \
	test_li_berger.cpp \
	test_linearalgebra.cpp \
	test_pohrt_li.cpp \
	test_main.cpp


OBJ = $(SRC:.cpp=.o)
TESTOBJ1 = $(TESTSRC:.cpp=.o)
TESTOBJ = $(TESTOBJ1:.cc=.o)

EXE = eval_stiffness eval_penetration eval_gap

unittests: $(TESTOBJ)
	$(LINK) $(LINKFLAGS) $(TESTOBJ) $(LIB) -o $@

$(EXE):	% : $(OBJ) %.o
	$(LINK) $(LINKFLAGS) $(OBJ) $@.o $(LIB) -o $@

%.o:%.cc
	$(CC) $(CCFLAGS) $(INC) -c $<

%.o:%.cpp
	$(CC) $(CCFLAGS) $(INC) -c $<

all: $(EXE)

clean:
	rm -f $(OBJ) $(TESTOBJ)
