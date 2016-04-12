/* ======================================================================
   USER-GFMD - Green's function molecular dynamics for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016) Lars Pastewka <lars.pastewka@kit.edu>, Tristan A. Sharp
   and others. See the AUTHORS file in the top-level USER-GFMD directory.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   ====================================================================== */
#ifndef __CRYSTAL_SURFACE_H
#define __CRYSTAL_SURFACE_H

#include <numeric>

#include "error.h"
#include "linearalgebra.h"

using namespace LAMMPS_NS;

/*!
 * Class describing surface unit cells, i.e. a unit cell that has one face
 * parallel to the surface. Its unit cell must be (a1,a2,a3) with a1=(100) and
 * a2=(010).
 */
class CrystalSurface {
 public:
    CrystalSurface();
    CrystalSurface(const CrystalSurface &);
    virtual ~CrystalSurface();

    /*!
     * Check whether the neighbors actually point to neighboring atoms
     */
    virtual bool check_neighbors(Error *error);

    struct lattice_neighbor_t {
        /*!
         * This atom's index within the neighbor structure.
         * Corresponds to the ordering in memory.
         */
        int indexj;

        /*!
         * Root atom's surface surface cell index.
         */
        int indexi;

        /*!
         * Type (in the LAMMPS sense) of this neighbor atom.
         */
        int jtype;

        /*!
         * Neighbors relative position [2pi/q] (since grid spacing = nn dist)
         */
        double r[3];
        
        // The following struct variables are derived from the above.

        /*!
         * Precomputed length of r
         */
        double rnorm;

        /*!
         * Neighbor shell, defined by distance and atom type.
         * Used to index the derivatives of the potential.
         * Convention: 0 for nn type, 1 for 2nn, 2 for 3nn
         * If two types of atoms are at the same distance, lower type
         * (of base atom) is listed first.
         * Could rederive this always from neighbor index.
         */
        int neighbor_shell;

        /*!
         * alpha minus beta (or "above"): rel. z-pos of this neighbor, units of
         * atomic planes
         */
        int x, y, ab;
    };

    const char *get_name() const {
        return name_;
    }

    const double *get_position(int i=0) const {
        return &positions_[3*i];
    }

    double *get_positions() const {
        return positions_;
    }

    /*!
     * The type of atom i; at which lattice neighbor struct is centered.
     */
    int get_type(int i=0) const {
        if (types_) {
            return types_[i];
        }
        return 1;
    }

    int *get_types() const {
        return types_;
    }

    void set_cell(const double *new_cell, Error *error) {
        // check whether first lattice vector is along x and second lattice
        // vector is along y
        if (new_cell[1] != 0.0 || new_cell[3] != 0.0 || new_cell[6] != 0.0 ||
            new_cell[7] != 0.0) {
            error->all(FLERR,"Cell must be aligned along x- and y-axis.");
        }
        memcpy(cell_, new_cell, 9*sizeof(double));
        compute_inverse_cell(error);
    }

    void set_cell_and_scale_atoms(const double *new_cell, Error *error);

    const double *get_cell() const {
        return cell_;
    }

    const double *get_invcell() const {
        return invcell_;
    }

    virtual int get_number_of_atoms() const {
        return num_atoms_;
    }

    virtual int get_number_of_neighbors(int i=-1) const {
        if (!lattice_neighbors_)
            return 0;
        if (i < 0) {
            return std::accumulate(num_neighbors_per_atom_,
                                   num_neighbors_per_atom_+num_atoms_,
                                   0);
        }
        else {
            return num_neighbors_per_atom_[i];
        }
    }

    virtual const lattice_neighbor_t *get_neighbors(int i=0) const {
        int k = 0;
        for (int j = 0; j < i; j++)  k += num_neighbors_per_atom_[j];
        return &lattice_neighbors_[k];
    }

    /*!
     * Get the distance between atom i, and the atom pointed at by n in terms
     * of cells that are spanned. Note: x can be fractional, in which case
     * something went wrong.
     */
    virtual void cell_distance(int i, const lattice_neighbor_t &n, double *x)
        const;

    /*!
     * Create a supercell of the lattice
     */
    virtual void supercell(int, int, int, Error *, int fx=-1, int fy=-1,
                           int fz=-1);

    /*!
     * Given a positions, compute the mapping of atoms in this positions to atoms in
     * the current positions.
     */
    virtual bool get_mapping(int, double *, int *, int *types=NULL,
                             double tol=1e-3) const;

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *f);

    /*!
     * Construct neighbor list up to a certain cutoff
     */
    void compute_neighbor_list(double);

    /*!
     * Reverse the order of atoms within the cell
     */
    void reverse_indices();

    /*!
     * Return cell index of cell within supercell
     */
    static int _cidx(int ny, int nz, int x, int y, int z) {
        return z+nz*(y+ny*x);
    }

    /*!
     * Return index of atom within supercell from index of atom within
     * reference cell
     */
    static int _idx(int ny, int nz, int ni, int x, int y, int z, int i,
                    int fx=-1, int fy=-1, int fz=-1) {
        int c = _cidx(ny, nz, x, y, z);
        int fc = _cidx(ny, nz, fx, fy, fz);

        // fc is cell whose atoms have smallest indices
        if (fc > 0) {
            if (c == fc)  c = 0;
            else if (c < fc)  c++;
        }

        return i+ni*c;
    }

 protected:
    void compute_inverse_cell(Error *);

    void identify_neighbor_shells();

    void transform(const double *);

    void _resize(int, int);

    /*!
     * Descriptor for this surface
     */
    char name_[80];

    /*!
     * Crystal unit cell and inverse (3, 3)
     */
    double cell_[9], invcell_[9];

    /*!
     * Number of atoms per surface cell.
     * The surface cell depth is the GF layer thickness.
     */
    int num_atoms_;

    /*!
     * Positions of the crystal (num_atoms_, 3)
     */
    double *positions_;

    /*!
     * Atom types as in the normal LAMMPS sense
     * Within a GF unit cell, can be multiple atom types
     * Atom types (num_atoms_, ) 
     */
    int *types_;

    /*!
     * Number of neighbors per positions atom
     */
    int *num_neighbors_per_atom_;

    /*!
     * Actual lattice neighbors
     */
    lattice_neighbor_t *lattice_neighbors_;
};

CrystalSurface *crystal_surface_factory(char *, int, int *, char **, Error *);


/*!
 * Overload operator<< for stream functions
 */
inline std::ostream &operator<<(std::ostream& os,
                                const CrystalSurface::lattice_neighbor_t &obj)
{
    os << "( " << obj.indexi << "-" << obj.indexj << ": r = (" << obj.r[0] 
       << "," << obj.r[1] << "," << obj.r[2] << "), |r| = " << obj.rnorm
       << ", x = (" << obj.x << "," << obj.y << "," << obj.ab << ") )";

    return os;
}

#endif
