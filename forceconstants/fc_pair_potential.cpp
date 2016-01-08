#include <algorithm>

#include <stdlib.h>
#include <string.h>

#include "linearalgebra.h"

#include "fc_pair_potential.h"

FCPairPotential::FCPairPotential(int narg, int *carg, char **arg,
                                 CrystalSurface *surface, Force *force,
                                 Error *error)
    : ForceConstants("pair-potential", NULL, error), linf_(NULL)
{
    const double *cell = surface->get_cell();
    surface->compute_neighbor_list(2.01*std::max(cell[0], cell[4]));

    // Input order: k1 (kappa1) k2 (kappa2) k3 (kappa3)
    k_[0] = 0.0;
    k_[1] = 0.0;
    k_[2] = 0.0;
    kappa_[0] = 0.0;
    kappa_[1] = 0.0;
    kappa_[2] = 0.0;

    if (*carg >= narg) {
        error->all(FLERR,"Missing argument: Number number of interacting "
                   "neighbors.");
    }

    // Before the the nu_ input location, find optional key word 'x' (ending in x)
    // Choosing to set kappa values
    char *s = strdup(arg[*carg]);

    bool x = false;
    int l = strlen(s);
    if (s[l-1] == 'x') {
        x = true;
        //printf("Setting the x option for kappa\n");
        s[l-1] = 0;
    }

    // Number of neighbors
    char *endptr;
    int nk = strtol(s, &endptr, 10); 
    //printf("nk is integer %s \n",arg[*carg]);
    if (endptr == arg[*carg]) {
        error->all(FLERR,"Error maximum neighbor shell to integer.");
    }
    (*carg)++;

    free(s);

    for (int i = 0; i < nk; i++) {
        if (*carg >= narg) {
            error->all(FLERR,"Spring constant expected.");
        }
        
        //printf("k_ is the float: %s \n",arg[*carg]);
        k_[i] = strtod(arg[*carg], &endptr);
        if (endptr == arg[*carg]) {
            error->all(FLERR,"Error converting spring constant to number.");
        }
        (*carg)++;
        
        if (x) {
            if (*carg >= narg) {
                error->all(FLERR,"Spring constant expected.");
            }
            
            //printf("kappa_ is the float: %s \n",arg[*carg]);
            kappa_[i] = strtod(arg[*carg], &endptr);
            if (endptr == arg[*carg]) {
                error->all(FLERR,"Error converting spring constant to "
                           "number.");
            }
            (*carg)++;
        }
    }

    bool keyword_found = true;
    while (*carg < narg-1 && keyword_found) {
        keyword_found = false;
        // add a linear force (per atom) contribution normal to the surface to
        // each layer to offset surface relaxation
        if (!strcmp(arg[*carg], "linf")) {
            // Number of layers
            (*carg)++;
            int nu = strtol(arg[*carg], &endptr, 10);
            if (endptr == arg[*carg]) {
                error->all(FLERR,"Error maximum number of layers to integer.");
            }

            linf_ = new double[nu];
            memset(linf_, 0, nu*sizeof(double));

            double linfsum = 0.0;
            for (int i = 0; i < nu; i++) {
                (*carg)++;
                if ((*carg) >= narg) {
                    char errstr[1024];
                    sprintf(errstr, "Stiffness kernel report %i layers,"
                            " please provide that many *linf* parameters.", nu);
                    error->all(FLERR, errstr);
                }
                linf_[i] = atof(arg[*carg]);
                linfsum += linf_[i];
            }
            if (abs(linfsum) > 1e-12) {
                char errstr[1024];
                sprintf(errstr, "Linear force contributions must sum to "
                        "zero but sums to %e.", linfsum);
                error->all(FLERR, errstr);
            }
            (*carg)++;
            keyword_found = true;
        }
    }
}


FCPairPotential::~FCPairPotential()
{
  if (linf_) {
    delete [] linf_;
  }
}


void FCPairPotential::offsite(const CrystalSurface *surf,
                              const CrystalSurface::lattice_neighbor_t *pair,
                              int ab, double *_D, Error *error)
{
    //    printf("FCPairPotential::get_force_constants\n");

    vec<double, 3> e(pair->r);
    e /= e.nrm2();

    mat<double> para(3, _D);
    outer(e, e, para);

    mat<double> perp = identity<double>(3) - para;

    int n = pair->neighbor_shell;
    double k = 0.0, kappa = 0.0;
    if (n < 3) {
        k = k_[n];
        kappa = kappa_[n];
    }

    //    printf("A\n");
    para *= -k;
    //    printf("B\n");
    perp *= kappa;
    //    printf("C\n");

    para -= perp;
    //    printf("D\n");
}


void FCPairPotential::linear(const CrystalSurface *surf, double *linf)
{
    if (linf_) {
        std::copy(linf_, linf_+surf->get_number_of_atoms(), linf);
    }
    else {
        std::fill(linf, linf+surf->get_number_of_atoms(), 0.0);
    }
}
