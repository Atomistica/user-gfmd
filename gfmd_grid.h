#ifndef GFMD_GRID_H
#define GFMD_GRID_H

/* ----------------------------------------------------------------------
 * Grid indices. The code uses a compound index based on a 32-bit or
 * 64-bit integer (depending on availability of that datatype).
 * --------------------------------------------------------------------*/

#define CONT_IDX(ix, iy, iu)  ( ( (ix)*(ny) + (iy) )*(nu) + (iu) )
#define LOC_IDX(ix, iy, iu)  ( ( (ix)*(ny_loc) + (iy) )*(nu) + (iu) )

#define IX_FROM_CONT_IDX(idx)  ( (idx) / ( (ny) * (nu) ) )
#define IY_FROM_CONT_IDX(idx)  ( ( (idx) / (nu) ) % (ny) )
#define IU_FROM_CONT_IDX(idx)  ( (idx) % (nu) )

#define POW2_TO_CONT_IDX(idx)  CONT_IDX(IX_FROM_POW2_IDX(idx), IY_FROM_POW2_IDX(idx), IU_FROM_POW2_IDX(idx))

#ifdef LAMMPS_SMALLSMALL

/*
 * bigint is a 32 bit integer:
 * This limits us to GFMD sizes of 8192 x 8192 and kernel sizes of 32
 * (10 atoms per unit cell) which should be okay for any conceivable
 * situation on a 32 bit machine.
 * IX, IY, IU are the grid indices.
 * IF is a flag, that if true, tells the pair style to ignore this atom
 * in the energy calculation.
 */
#define POW2_IDX(ix, iy, iu, flag) ( (bigint(ix) << 19) + (bigint(iy) << 6) + (bigint(iu) << 1) + flag )

#define IX_FROM_POW2_IDX(idx)  ( idx >> 19 )
#define IY_FROM_POW2_IDX(idx)  ( ( idx >> 6 ) & 8191 )
#define IU_FROM_POW2_IDX(idx)  ( ( idx >> 1 ) & 31 )
#define FLAG_FROM_POW2_IDX(idx)  ( idx & 1 )

#else

/*
 * bigint is a 64 bit integer:
 * This limits us to GFMD sizes of 67108864 x 67108864 and kernel sizes of 2048
 * (682 atoms per unit cell) which should be okay for any conceivable
 * situation on a 64 bit machine
 * IX, IY, IU are the grid indices.
 * FLAG, if true, tells the pair style to ignore this atom
 * in the energy calculation.
 */
#define POW2_IDX(ix, iy, iu, flag) ( (bigint(ix) << 38) + (bigint(iy) << 12) + (bigint(iu) << 1) + flag )

#define IX_FROM_POW2_IDX(idx)  ( idx >> 38 )
#define IY_FROM_POW2_IDX(idx)  ( ( idx >> 12 ) & 67108863 )
#define IU_FROM_POW2_IDX(idx)  ( ( idx >> 1 ) & 2047 )
#define FLAG_FROM_POW2_IDX(idx)  ( idx & 1 )

#endif

#endif
