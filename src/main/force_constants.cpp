/* ======================================================================
   USER-GFMD - Elastic half-space methods for LAMMPS
   https://github.com/Atomistica/user-gfmd

   Copyright (2011-2016,2021)
      Lars Pastewka <lars.pastewka@imtek.uni-freiburg>,
      Tristan A. Sharp and others.
   See the AUTHORS file in the top-level USER-GFMD directory.

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
#include <string.h>

#include "force_constants.h"
#ifndef PYTHON
#include "fc_lj_cut.h"
#include "fc_lj_cut_fd.h"
#include "fc_lj_smooth.h"
#include "fc_pair_potential.h"
#endif

#ifdef GFMD_EAM
#include "fc_eam.h"
#include "fc_eam_fd.h"
#endif

#ifdef GFMD_TERSOFF
#include "fc_tersoff.h"
#endif

#ifdef GFMD_ATOMISTICA
#include "fc_atomistica.h"
#endif

ForceConstants::ForceConstants(const char *name, Force *force, Error *error,
                               const char *pair_style) : manybody_(false)
{
    const char *ps = pair_style;
    strcpy(name_, name);
    if (!ps) {
        ps = name_;
    }

    if (force) {
        if (strcmp(force->pair_style, ps)) {
            char errstr[1024];
            sprintf(errstr, "Declared force constants pair style '%s' does "
                    "not agree with pair style '%s'.", ps, force->pair_style);
            error->all(FLERR,errstr);
        }
    }
}


/*!
 * Dump generic text info to file
 */
void ForceConstants::dump_info(FILE *f, CrystalSurface *surf, Error *error)
{
    fprintf(f, "Interaction model = '%s'.\n", name_);
}


/*!
 * On-site contribution to the force constant matrix, i.e. the
 * 3*n_atoms x 3*n_atoms matrix D_ij for i == j.
 */
void ForceConstants::onsite(const CrystalSurface *surf, int i, int ab,
                            double *_D, Error *error)
{
    mat<double> D(3, _D);
    D = 0.0;

    int n = surf->get_number_of_neighbors(i);
    const CrystalSurface::lattice_neighbor_t *neigh = surf->get_neighbors(i);

    for (int j = 0; j < n; j++, neigh++) {
        // For the topmost layer, we need to exclude parts of the sum.
        // These are the U0, U1, etc. surface matrices.
        if (neigh->ab <= ab) {
            double tmp[9];
            offsite(surf, neigh, ab, tmp, error);
            D -= tmp;
        }
    }
}


/*!
 * Linear force contribution.
 */
void ForceConstants::linear(const CrystalSurface *surf, double *linf,
                            Error *error)
{
    std::fill(linf, linf+surf->get_number_of_atoms(), 0.0);
}


/*!
 * Instantiate a ForceConstants object according to keyword arguments
 */
ForceConstants *force_constants_factory(char *keyword, int narg, int *carg,
                                        char **arg, CrystalSurface *surface,
                                        Force *force, Memory *memory,
                                        Error *error)
{
    char errstr[1024];
    ForceConstants *fc = NULL;

#define FC_CLASS
#define FCStyle(key, Class)                                              \
    if (!strcmp(keyword, #key)) {                                        \
        fc = new Class(narg, carg, arg, surface, force, error);         \
    }
#define FCStyleWithMemory(key, Class)                                   \
    if (!strcmp(keyword, #key)) {                                        \
        fc = new Class(narg, carg, arg, surface, force, memory, error);  \
    }

#ifndef PYTHON
#include "fc_lj_cut.h"
#include "fc_lj_cut_fd.h"
#include "fc_lj_smooth.h"
#include "fc_pair_potential.h"
#endif

#ifdef GFMD_EAM
#include "fc_eam.h"
#include "fc_eam_fd.h"
#endif

#ifdef GFMD_TERSOFF
#include "fc_tersoff.h"
#endif

#ifdef GFMD_ATOMISTICA
#include "fc_atomistica.h"
#endif

#undef FCStyle
#undef FCStyleWithMemory

    if (fc) {
        if (strcmp(fc->get_name(), keyword)) {
            sprintf(errstr, "force_constants_factory: Internal error: "
                    "keyword '%s' and ForceConstants '%s' name mismatch.",
                    keyword, fc->get_name());
            error->all(FLERR,errstr);
        }
    }

    return fc;
}
