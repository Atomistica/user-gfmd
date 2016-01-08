CXXFLAGS = -g -O3 -funroll-loops

bicubic.so:  memory.o linearalgebra.o table2d.o py_bicubic.o
	g++ ${CXXFLAGS} -shared memory.o linearalgebra.o table2d.o py_bicubic.o -o bicubic.so


surface_stiffness.so: memory.o linearalgebra.o surface_stiffness.o debug_stiffness.o isotropic_stiffness.o sc100_stiffness.o fcc100_stiffness.o py_surface_stiffness.o
	g++ ${CXXFLAGS} -shared memory.o linearalgebra.o surface_stiffness.o debug_stiffness.o isotropic_stiffness.o sc100_stiffness.o fcc100_stiffness.o py_surface_stiffness.o -o surface_stiffness.so


py_bicubic.o:  py_bicubic.cpp py_bicubic.h
	g++  ${CXXFLAGS} -fPIC -c py_bicubic.cpp `python-config --cflags` -I`python -c "import numpy ; print numpy.get_include()"` -I.. -I../LAMMPS_STUB


py_surface_stiffness.o:  py_surface_stiffness.cpp py_surface_stiffness.h
	g++  ${CXXFLAGS} -fPIC -c py_surface_stiffness.cpp `python-config --cflags` -I`python -c "import numpy ; print numpy.get_include()"` -I.. -I../LAMMPS_STUB


memory.o:  ../LAMMPS_STUB/memory.cpp ../LAMMPS_STUB/memory.h
	g++  ${CXXFLAGS} -fPIC -I../LAMMPS_STUB -c ../LAMMPS_STUB/memory.cpp

linearalgebra.o:  ../linearalgebra.cpp ../linearalgebra.h
	g++  ${CXXFLAGS} -fPIC -I../LAMMPS_STUB -c ../linearalgebra.cpp

table2d.o:  ../table2d.cpp ../table2d.h
	g++  ${CXXFLAGS} -fPIC -I../LAMMPS_STUB -c ../table2d.cpp

surface_stiffness.o:  ../surface_stiffness.cpp ../surface_stiffness.h
	g++  ${CXXFLAGS} -fPIC -I../LAMMPS_STUB -c ../surface_stiffness.cpp

debug_stiffness.o:  ../debug_stiffness.cpp ../debug_stiffness.h
	g++  ${CXXFLAGS} -fPIC -I../LAMMPS_STUB -c ../debug_stiffness.cpp

isotropic_stiffness.o:  ../isotropic_stiffness.cpp ../isotropic_stiffness.h
	g++  ${CXXFLAGS} -fPIC -I../LAMMPS_STUB -c ../isotropic_stiffness.cpp

sc100_stiffness.o:  ../sc100_stiffness.cpp ../sc100_stiffness.h
	g++  ${CXXFLAGS} -fPIC -I../LAMMPS_STUB -c ../sc100_stiffness.cpp

fcc100_stiffness.o:  ../fcc100_stiffness.cpp ../fcc100_stiffness.h
	g++  ${CXXFLAGS} -fPIC -I../LAMMPS_STUB -c ../fcc100_stiffness.cpp


all:	bicubic.so surface_stiffness.so


clean:
	rm *.so *.o
