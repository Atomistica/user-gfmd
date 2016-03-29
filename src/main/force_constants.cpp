#include <string.h>

#include "force_constants.h"
#include "fc_lj_cut.h"
#include "fc_lj_cut_fd.h"
#include "fc_lj_smooth.h"
#include "fc_pair_potential.h"

#ifdef GFMD_MANYBODY
#include "fc_eam.h"
#include "fc_eam_fd.h"
#include "fc_tersoff.h"
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
    char errstr[120];
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

#include "fc_lj_cut.h"
#include "fc_lj_cut_fd.h"
#include "fc_lj_smooth.h"
#include "fc_pair_potential.h"

#ifdef GFMD_MANYBODY
#include "fc_eam.h"
#include "fc_eam_fd.h"
#include "fc_tersoff.h"
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
