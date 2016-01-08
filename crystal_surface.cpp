#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

#include "linearalgebra.h"

#include "crystal_surface.h"
#include "dia100_surface.h"
#include "dia111_surface.h"
#include "fcc100_surface.h"
#include "fcc111_surface.h"
#include "sc100_surface.h"


CrystalSurface::CrystalSurface()
    : positions_(NULL), types_(NULL), lattice_neighbors_(NULL),
      num_neighbors_per_atom_(NULL)
{
}


CrystalSurface::CrystalSurface(const CrystalSurface &other)
    : positions_(NULL), types_(NULL), lattice_neighbors_(NULL),
      num_neighbors_per_atom_(NULL)
{
    strcpy(name_, other.name_);
    std::copy(other.cell_, other.cell_+9, cell_);
    std::copy(other.invcell_, other.invcell_+9, invcell_);
    num_atoms_ = other.num_atoms_;
    positions_ = new double[3*num_atoms_];
    std::copy(other.positions_, other.positions_+3*num_atoms_, positions_);
    types_ = new int[num_atoms_];
    std::copy(other.types_, other.types_+num_atoms_, types_);
    int num_neighb = other.get_number_of_neighbors();
    if (num_neighb > 0) {
        num_neighbors_per_atom_ = new int[num_atoms_];
        std::copy(other.num_neighbors_per_atom_,
                  other.num_neighbors_per_atom_+num_atoms_,
                  num_neighbors_per_atom_);   

        lattice_neighbors_ = new lattice_neighbor_t[num_neighb];
        std::copy(other.lattice_neighbors_, other.lattice_neighbors_+num_neighb,
                  lattice_neighbors_);
        assert(get_number_of_neighbors() == num_neighb);
    }
}


CrystalSurface::~CrystalSurface()
{
    if (positions_) {
        delete [] positions_;
    }

    if (types_) {
        delete [] types_;
    }

    if (num_neighbors_per_atom_) {
        delete [] num_neighbors_per_atom_;
    }

    if (lattice_neighbors_) {
        delete [] lattice_neighbors_;
    }
}


/*!
 * Dump generic text info to file
 */
void CrystalSurface::dump_info(FILE *f)
{
    fprintf(f, "Surface = '%s'. Surface unit cell definition follows.\n",
            name_);
    fprintf(f, "Lattice vectors:\n");
    fprintf(f, "  / %f \\   / %f \\   / %f \\\n", cell_[0], cell_[1],
            cell_[2]);
    fprintf(f, "  | %f |,  | %f |,  | %f |\n", cell_[3], cell_[4], cell_[5]);
    fprintf(f, "  \\ %f /   \\ %f /   \\ %f /\n", cell_[6], cell_[7],
            cell_[8]);
    fprintf(f, "Atomic positions [fractional atomic positions]:\n");
    for (int i = 0; i < num_atoms_; i++) {
        double x[3];
        mat_mul_vec(3, invcell_, &positions_[3*i], x);
        fprintf(f, "  u=%i --- ( %f, %f, %f )  [ ( %f, %f, %f ) ]\n", i,
                positions_[3*i+0], positions_[3*i+1], positions_[3*i+2], x[0], 
                x[1], x[2]);
    }

    // Count number of neighbors per shell.
    const int shellmax = 3;
    int neighbors_per_shell[shellmax*num_atoms_];
    std::fill(neighbors_per_shell, neighbors_per_shell+shellmax*num_atoms_, 0);

    int nn = get_number_of_neighbors();
    const CrystalSurface::lattice_neighbor_t *neigh = lattice_neighbors_;
    for (int i = 0; i < nn; i++, neigh++) {
        int ns = neigh->neighbor_shell;
        if (ns < shellmax) {
            neighbors_per_shell[ns*num_atoms_+neigh->indexi]++;
        }
    }

    // Out number of neighbors per shell.
    fprintf(f, "Number of neighbors (neighbor shell):\n");
    for (int i = 0; i < num_atoms_; i++) {
        fprintf(f, "  u=%i ---", i);
        for (int j = 0; j < shellmax; j++) {
            if (neighbors_per_shell[j*num_atoms_+i] > 0) {
                fprintf(f, " %i(%i)", neighbors_per_shell[j*num_atoms_+i],
                        j+1);
            }
        }
        fprintf(f, "\n");
    }
}


void CrystalSurface::compute_inverse_cell(Error *error)
{
    //conj_transpose(3, invcell_, cell_);
    memcpy(invcell_, cell_, 9*sizeof(double));
    GaussJordan(3, invcell_, error);
}


void CrystalSurface::_resize(int old_size, int new_size)
{
    lattice_neighbor_t *new_lattice_neighbors = 
        new lattice_neighbor_t[new_size];
    if (lattice_neighbors_) {
        std::copy(lattice_neighbors_,
                  lattice_neighbors_+std::min(old_size, new_size),
                  new_lattice_neighbors);
        delete [] lattice_neighbors_;
    }
    lattice_neighbors_ = new_lattice_neighbors;
}


void CrystalSurface::compute_neighbor_list(double rcut)
{
    if (num_neighbors_per_atom_)
        delete [] num_neighbors_per_atom_;
    num_neighbors_per_atom_ = new int[num_atoms_];

    double lx = cell_[0]*cell_[0] + cell_[3]*cell_[3] + cell_[6]*cell_[6];
    double ly = cell_[1]*cell_[1] + cell_[4]*cell_[4] + cell_[7]*cell_[7];
    double lz = cell_[2]*cell_[2] + cell_[5]*cell_[5] + cell_[8]*cell_[8];

    int nx = lrint(floor(rcut/lx))+1;
    int ny = lrint(floor(rcut/ly))+1;
    int nz = lrint(floor(rcut/lz))+1;

    int nsize = 100;
    _resize(0, nsize);

    mat<double,3> cell(cell_);

    int n = 0;
    for (int i = 0; i < num_atoms_; i++) {
        int nn = 0;
        vec<double> ri(3, &positions_[3*i]);
        for (int j = 0; j < num_atoms_; j++) {
            for (int x = -nx; x <= nx; x++) {
                for (int y = -ny; y <= ny; y++) {
                    for (int z = -nz; z <= nz; z++) {
                        if (i != j || x != 0 || y != 0 || z != 0) {
                            vec<double> r(3);
                            r = ri;
                            r -= &positions_[3*j];

                            double d[3] = { x, y, z };
                            r -= cell.dot(d);

                            double rnorm = r.nrm2();

                            if (rnorm < rcut) {
                                if (n >= nsize) {
                                    _resize(nsize, 2*nsize);
                                    nsize *= 2;
                                }

                                lattice_neighbor_t *neigh = &lattice_neighbors_[n];
                                neigh->indexi = i;
                                neigh->indexj = j;
                                
                                vec<double>(3, neigh->r) = r;
                                neigh->rnorm = rnorm;

                                neigh->x = x;
                                neigh->y = y;
                                neigh->ab = z;

                                n++;
                                nn++;
                            }
                        }
                    }
                }
            }
        }

        num_neighbors_per_atom_[i] = nn;
    }

    identify_neighbor_shells();
}


void CrystalSurface::identify_neighbor_shells()
{
    int nr = 0;
    std::vector<double> r(1024);
    int c[1024];

    int num_neighbors = std::accumulate(num_neighbors_per_atom_,
                                        num_neighbors_per_atom_+num_atoms_,
                                        0);

    lattice_neighbor_t *neigh = lattice_neighbors_;
    for (int i = 0; i < num_neighbors; i++) {
        bool found = false;
        for (int k = 0; k < nr && !found; k++) {
            if (std::abs(neigh->rnorm-r[k]) < 0.01) {
                c[k]++;
                found = true;
            }
        }
        if (!found) {
            r[nr] = neigh->rnorm;
            c[nr] = 1;
            nr++;
        }
        neigh++;
    }

    r.resize(nr);
    std::sort(r.begin(), r.end());

    neigh = lattice_neighbors_;
    for (int i = 0; i < num_neighbors; i++) {
        bool found = false;
        for (int k = 0; k < nr && !found; k++) {
            if (std::abs(neigh->rnorm-r[k]) < 0.01) {
                neigh->neighbor_shell = k;
                found = true;
            }
        }

        neigh++;
    }
}


void CrystalSurface::reverse_indices()
{
    // New descriptors.

    double *new_positions = new double[3*num_atoms_];
    int *new_types = new int[num_atoms_];

    int *new_num_neighbors_per_atom = NULL;
    int num_lattice_neighbors = 0;
    if (num_neighbors_per_atom_) {
        new_num_neighbors_per_atom = new int[num_atoms_];
        num_lattice_neighbors = 
            std::accumulate(num_neighbors_per_atom_,
                            num_neighbors_per_atom_+num_atoms_,
                            0);
    }

    lattice_neighbor_t *new_lattice_neighbors = NULL;

    if (lattice_neighbors_) {
        new_lattice_neighbors = 
            new lattice_neighbor_t[num_lattice_neighbors];
        memset(new_lattice_neighbors, 0,
               num_lattice_neighbors*sizeof(lattice_neighbor_t));
    }

    // Compute new positions and copy lattice neighbors.
    lattice_neighbor_t *new_neigh = new_lattice_neighbors;
    for (int i = 0; i < num_atoms_; i++) {
        int k = num_atoms_-1-i;
        // Copy positions.
        new_positions[3*i] = positions_[3*k];
        new_positions[3*i+1] = positions_[3*k+1];
        new_positions[3*i+2] = positions_[3*k+2];

        // Copy type.
        new_types[i] = types_[k];

        // Copy lattice neighbors.
        if (num_neighbors_per_atom_) {
            new_num_neighbors_per_atom[i] = num_neighbors_per_atom_[k];
            if (lattice_neighbors_) {
                const lattice_neighbor_t *neigh = get_neighbors(k);
                for (int n = 0; n < num_neighbors_per_atom_[k]; n++) {
                    new_neigh->r[0] = neigh->r[0];
                    new_neigh->r[1] = neigh->r[1];
                    new_neigh->r[2] = neigh->r[2];
                    new_neigh->rnorm = neigh->rnorm;
                    new_neigh->neighbor_shell = neigh->neighbor_shell;
                    
                    new_neigh->indexi = num_atoms_-1-neigh->indexi;
                    new_neigh->indexj = num_atoms_-1-neigh->indexj;

                    new_neigh->x = neigh->x;
                    new_neigh->y = neigh->y;
                    new_neigh->ab = neigh->ab;
                    
                    neigh++;
                    new_neigh++;
                }
            }
        }
    }

    // Delete old, and copy new to old.

    delete [] positions_;
    delete [] types_;
    if (num_neighbors_per_atom_)
        delete [] num_neighbors_per_atom_;
    if (lattice_neighbors_)
        delete [] lattice_neighbors_;

    positions_ = new_positions;
    types_ = new_types;
    if (num_neighbors_per_atom_)
        num_neighbors_per_atom_ = new_num_neighbors_per_atom;
    if (lattice_neighbors_)
        lattice_neighbors_ = new_lattice_neighbors;
}


void CrystalSurface::transform(const double *_T)
{
    mat<double> T(3, _T, true);

    for (int i = 0; i < num_atoms_; i++) {
        vec<double> position(3, &positions_[3*i]);
        position = T.dot(position);
    }

    if (lattice_neighbors_) {
        lattice_neighbor_t *neigh = lattice_neighbors_;
        for (int i = 0; i < num_atoms_; i++) {
            for (int n = 0; n < num_neighbors_per_atom_[i]; n++) {
                vec<double> r(3, neigh->r);
                r = T.dot(r);
                neigh->rnorm = r.nrm2();
                neigh++;
            }
        }   
    }
}


void CrystalSurface::set_cell_and_scale_atoms(const double *new_cell, Error *error) {
    // check whether first lattice vector is along x and second lattice
    // vector is along y
    if (new_cell[1] != 0.0 || new_cell[3] != 0.0 || new_cell[6] != 0.0 ||
        new_cell[7] != 0.0) {
        error->all(FLERR,"Cell must be aligned along x- and y-axis.");
    }

    // Convert all vectors to scaled coordinates
    transform(invcell_);

    // Copy new cell information and compute inverse
    memcpy(cell_, new_cell, 9*sizeof(double));
    compute_inverse_cell(error);

    // Convert all vectors to unscaled coordinates
    transform(cell_);
}


void CrystalSurface::cell_distance(int i, const lattice_neighbor_t &n,
                                   double *x) const
{
    // return r_i + n.r - r_j, scaled by the inverse lattice
    assert(i == n.indexi);

    // compute distance between atoms in cell
    vec<double,3> r = &positions_[3*i];
    r -= &positions_[3*n.indexj];

    // compute distance to neighbor minus any offset within the cell
    r -= n.r;

    // x is the scaled lattice coordinate
    mat_mul_vec(3, invcell_, r.data(), x);
}


void CrystalSurface::supercell(int nx, int ny, int nz, Error *error, int fx,
                               int fy, int fz)
{
    // new descriptors

    assert(nx > 0);
    assert(ny > 0);
    assert(nz > 0);

    int num_lattice_neighbors = get_number_of_neighbors();

    int new_num_atoms = num_atoms_*nx*ny*nz;

    double *new_positions = new double[3*new_num_atoms];
    int *new_types = new int[new_num_atoms];

    int *new_num_neighbors_per_atom = NULL;
    int new_num_lattice_neighbors = 0;
    if (num_neighbors_per_atom_) {
        new_num_neighbors_per_atom = new int[new_num_atoms];
        new_num_lattice_neighbors = nx*ny*nz*num_lattice_neighbors;
    }

    lattice_neighbor_t *new_lattice_neighbors = NULL;

    if (lattice_neighbors_) {
        new_lattice_neighbors = 
            new lattice_neighbor_t[new_num_lattice_neighbors];
    }

    // Compute new cell.

    mat<double,3> cell(cell_);
    mat<double,3> new_cell(cell_);

    new_cell[0][0] *= nx;
    new_cell[1][0] *= nx;
    new_cell[2][0] *= nx;
    new_cell[0][1] *= ny;
    new_cell[1][1] *= ny;
    new_cell[2][1] *= ny;
    new_cell[0][2] *= nz;
    new_cell[1][2] *= nz;
    new_cell[2][2] *= nz;

    // Make sure that new cell has minimal possible shift in x- and y-
    // direction.
    new_cell[0][2] -= floor(new_cell[0][2]/new_cell[0][0])*new_cell[0][0];
    new_cell[1][2] -= floor(new_cell[1][2]/new_cell[1][1])*new_cell[1][1];

    // Set cell and compute inverse cell.
    set_cell(new_cell.const_data(), error);

    // Compute new positions and copy lattice neighbors.
    for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
    for (int z = 0; z < nz; z++) {
        int k = _idx(ny, nz, num_atoms_, x, y, z, 0, fx, fy, fz);

        lattice_neighbor_t *neigh = lattice_neighbors_;
        lattice_neighbor_t *new_neigh = new_lattice_neighbors;
        if (new_neigh) {
            new_neigh += k*num_lattice_neighbors/num_atoms_;
        }

        for (int i = 0; i < num_atoms_; i++, k++) {
            // wrap positions atom into cell
            vec<double,3> pos(&positions_[3*i]);
            vec<double,3> c;
            c[0] = x;
            c[1] = y;
            c[2] = z;
            pos += cell.dot(c);
            double spos[3];
            mat_mul_vec(3, invcell_, pos.const_data(), spos);
            spos[0] -= floor(spos[0]);
            spos[1] -= floor(spos[1]);
            spos[2] -= floor(spos[2]);
            mat_mul_vec(3, cell_, spos, &new_positions[3*k]);

            // set type
            new_types[k] = types_[i];

            // copy lattice neighbors
            if (neigh) {
                new_num_neighbors_per_atom[k] = num_neighbors_per_atom_[i];
                for (int n = 0; n < num_neighbors_per_atom_[i]; n++) {
                    new_neigh->r[0] = neigh->r[0];
                    new_neigh->r[1] = neigh->r[1];
                    new_neigh->r[2] = neigh->r[2];
                    new_neigh->rnorm = neigh->rnorm;
                    new_neigh->neighbor_shell = neigh->neighbor_shell;

                    assert(neigh->indexi == i);
                    assert(k == _idx(ny, nz, num_atoms_, x, y, z, i, fx, fy,
                                     fz));
                    new_neigh->indexi = k;

                    // Compound_index is a continuous atoms index, starting at
                    // atom 0 of cell with nu = 0.
                    int new_index = _idx(ny, nz, num_atoms_,
                                         modulo(x+neigh->x, nx),
                                         modulo(y+neigh->y, ny),
                                         modulo(z+neigh->ab, nz),
                                         neigh->indexj,
                                         fx, fy, fz);

                    // new_index is index within new supercell.
                    new_neigh->indexj = new_index;

                    // This is new lattice offset,
                    // i.e. compound_index/new_num_atoms.
                    new_neigh->x = 
                        int(nearbyint(std::floor(double(x+neigh->x)/nx)));
                    new_neigh->y = 
                        int(nearbyint(std::floor(double(y+neigh->y)/ny)));
                    new_neigh->ab = 
                        int(nearbyint(std::floor(double(z+neigh->ab)/nz)));

                    neigh++;
                    new_neigh++;
                }
            }
        }
    }
    }
    }

    // Delete old, and copy new to old.

    delete [] positions_;
    delete [] types_;
    if (num_neighbors_per_atom_)
        delete [] num_neighbors_per_atom_;
    if (lattice_neighbors_)
        delete [] lattice_neighbors_;

    num_atoms_ = new_num_atoms;
    positions_ = new_positions;
    types_ = new_types;
    if (new_num_neighbors_per_atom)
        num_neighbors_per_atom_ = new_num_neighbors_per_atom;
    if (new_lattice_neighbors)
        lattice_neighbors_ = new_lattice_neighbors;
}


bool CrystalSurface::get_mapping(int num_ref, double *ref_positions,
                                 int *indices, int *types, double tol) const
{
    bool mapping_found = false;

    // i is our guess for which atom corresponds to the first unit cell atom.
    for (int i = 0; i < num_atoms_ && !mapping_found; i++) {

        // Initialize mapping.
        for (int j = 0; j < num_ref; j++) {
            indices[j] = -1;
        }

        indices[i] = 0;
        types_[0] = types[i];
        
        // check distance to all other atoms
        for (int j = 0; j < num_ref; j++) {
            if (i != j) {
                // dr1 is distance between atom i and j
                vec<double,3> dr1(&ref_positions[3*i]);
                dr1 -= &ref_positions[3*j];

                // Now check which distance this corresponds to.
                for (int k = 0; k < num_atoms_; k++) {
                    vec<double,3> dr2(&positions_[0]);
                    dr2 -= &positions_[3*k];

                    // Are distance vectors identical?
                    if (dr1.almost_equal(dr2, tol)) {
                        indices[j] = k;
                        if (types) {
                            types_[k] = types[j];
                        }
                    }
                }
            }
        }

        // Check if we have num_atoms_ indices.
        int num_found = num_ref - std::count(indices, indices+num_ref, -1);
        mapping_found = num_found == num_atoms_;
    }

    if (!mapping_found) {
        std::cout << "Unit cell positions:" << std::endl;
        for (int i = 0; i < num_atoms_; i++) {
            std::cout << vec<double>(3, &positions_[3*i]) << std::endl;
        }
        std::cout << "Reference positions:" << std::endl;
        for (int i = 0; i < num_ref; i++) {
            std::cout << vec<double>(3, &ref_positions[3*i]) << std::endl;
        }
    }

    return mapping_found;
}


/*!
 * Check whether the neighbors actually point to neighboring atoms
 */
bool CrystalSurface::check_neighbors(Error *error)
{
    bool res = true;

    int k = 0;
    for (int i = 0; i < num_atoms_; i++, k += num_neighbors_per_atom_[i]) {
        for (int j = 0; j < num_neighbors_per_atom_[i]; j++) {
            // x is the scaled lattice coordinate
            double x[3];
            cell_distance(i, lattice_neighbors_[k+j], x);

            // x should be three integer numbers
            for (int l = 0; l < 3; l++) {
                if (std::abs(x[l]-round(x[l])) > 1e-6) {
                    std::cout << "l = " << l << ", x[l]-round(x[l]) = " << 
                        x[l]-round(x[l]) << ", x[l] = " << x[l] << std::endl;
                    res = false;
                }
            }

            // should be x distance (in number of cells)
            if (std::abs(round(x[0])-lattice_neighbors_[k+j].x) > 1e-6) {
                std::cout << "j = " << j << ", k = " << k << 
                    ", round(x[0])-lattice_neighbors_[k+j].x = " << 
                    round(x[0])-lattice_neighbors_[k+j].x << ",  x[0] = " <<
                    x[0] << ", lattice_neighbors_[k+j].x = " <<
                    lattice_neighbors_[k+j].x << std::endl;
                res = false;
            }

            // should be y distance (in number of cells)
            if (std::abs(round(x[1])-lattice_neighbors_[k+j].y) > 1e-6) {
                std::cout << "j = " << j << ", k = " << k << 
                    ", round(x[1])-lattice_neighbors_[k+j].y = " << 
                    round(x[1])-lattice_neighbors_[k+j].y << ",  x[1] = " <<
                    x[1] << ", lattice_neighbors_[k+j].y = " <<
                    lattice_neighbors_[k+j].y << std::endl;
                res = false;
            }

            // ab should be z distance (in number of cells)
            if (std::abs(round(x[2])-lattice_neighbors_[k+j].ab) > 1e-6) {
                std::cout << "j = " << j << ", k = " << k << 
                    ", round(x[2])-lattice_neighbors_[k+j].ab = " << 
                    round(x[2])-lattice_neighbors_[k+j].ab << ",  x[2] = " <<
                    x[2] << ", lattice_neighbors_[k+j].ab = " <<
                    lattice_neighbors_[k+j].ab << std::endl;
                res = false;
            }
        }
    }

    return res;
}


/*!
 * Instantiate a CrystalSurface object according to keyword arguments
 */
CrystalSurface *crystal_surface_factory(char *keyword, int narg, int *carg,
                                        char **arg, Error *error)
{
    char errstr[120];
    CrystalSurface *surface = NULL;

#define CRYSTAL_SURFACE_CLASS
#define CrystalSurfaceStyle(key, Class)              \
    if (!strcmp(keyword, #key)) {                    \
        surface = new Class(narg, carg, arg, error); \
    }
#include "dia100_surface.h"
#include "dia111_surface.h"
#include "fcc100_surface.h"
#include "fcc111_surface.h"
#include "sc100_surface.h"
#undef CrystalSurfaceStyle

  if (surface) {
    if (strcmp(surface->get_name(), keyword)) {
      sprintf(errstr, "crystal_surface_factory: Internal error: keyword '%s' "
              "and CrystalSurface '%s' name mismatch.", keyword,
              surface->get_name());
      error->all(FLERR,errstr);
    }
  }

  return surface;
}

