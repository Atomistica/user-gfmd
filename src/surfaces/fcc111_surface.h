#ifdef CRYSTAL_SURFACE_CLASS

CrystalSurfaceStyle(fcc111,FCC111Surface)

#else

#ifndef __FCC111_SURFACE_H
#define __FCC111_SURFACE_H

#include "crystal_surface.h"

class FCC111Surface : public CrystalSurface {
 public:
	FCC111Surface(int, int *, char **, Error *);

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
