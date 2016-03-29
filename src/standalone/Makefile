SHELL      = /bin/sh

VPATH      = ../LAMMPS_STUB
VPATH     += ..

DEFINES    = -DNO_LAMMPS
CC         = icpc
CCFLAGS    = -O3 -openmp `nc-config --cflags` -I/home/pas/local/sdl-1.2.14/include/SDL
#CCFLAGS    = -O3 `nc-config --cflags` -I/home/pas/local/sdl-1.2.14/include/SDL
#CCFLAGS    = -O0 -g -pg `nc-config --cflags`
LINK       = icpc
LINKFLAGS  = -O3 -openmp
#LINKFLAGS  = -O3
#LINKFLAGS  = -O0 -g -pg

INC        = -I../LAMMPS_STUB -I.. -I.
LIB        =  `nc-config --libs` -L/home/pas/local/sdl-1.2.14/lib -lSDLmain -lSDL -lnetcdf_c++ -L/home/pas/local/fftw-3.3/lib -lfftw3_omp -lfftw3
#LIB        =  `nc-config --libs` -L/home/pas/local/sdl-1.2.14/lib -lSDLmain -lSDL -lnetcdf_c++ -lfftw3

SRC        = \
	mpi.cpp \
	memory.cpp \
	linearalgebra.cpp \
	atom.cpp \
	gfmd_memory.cpp \
	gfmd_stiffness.cpp \
	gfmd_misc.cpp \
	table2d.cpp \
	fix_contact_map.cpp \
	nc_traj_io.cpp

OBJ        = $(SRC:.cpp=.o)

EXE        = contact

$(EXE):	% : $(OBJ) %.o
	$(LINK) $(OBJ) $@.o $(LINKFLAGS) $(LIB) -o $@

%.o:%.cpp
	$(CC) $(CCFLAGS) $(DEFINES) $(INC) -c $< -o $@

all: $(EXE)

clean:
	rm $(OBJ) contact.o
