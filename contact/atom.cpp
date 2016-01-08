#include "atom.h"

using namespace LAMMPS_NS;

Atom::Atom(bigint in_natoms, Memory *in_memory)
{
  memory = in_memory;

  natoms = in_natoms;
  nlocal = natoms;
  nghost = 0;
  nmax = natoms;

  memory->create(x, nmax, 3, "Atom::x");
  memory->create(v, nmax, 3, "Atom::v");
  memory->create(f, nmax, 3, "Atom::f");
  memory->create(xeq, nmax, 3, "Atom::xeq");
  memory->create(mask, nmax, "Atom::mask");

  for (int i = 0; i < nmax; i++) {
    mask[i] = 1;
  }
}


Atom::~Atom()
{
  memory->destroy(x);
  memory->destroy(v);
  memory->destroy(f);
  memory->destroy(xeq);
}
