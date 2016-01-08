#ifdef CRYSTAL_SURFACE_CLASS

CrystalSurfaceStyle(dia111,Diamond111Surface)

#else

#ifndef __DIAMOND111_SURFACE_H
#define __DIAMOND111_SURFACE_H

#include "crystal_surface.h"

class Diamond111Surface : public CrystalSurface {
 public:
	Diamond111Surface(int, int *, char **, Error *);
	Diamond111Surface(double, int, Error *);

    /*!
     * Dump generic text info to file
     */
    virtual void dump_info(FILE *f);

 protected:
    /*!
     * Create primitive unit cell
     */
    void primitive_cell(Error *error);

    /*!
     * Number of unit cells (along surface normal).
     */
    int nu_;

    /*!
     * The lattice constant.
     */
    double lattice_constant_;
};

#endif

#endif
