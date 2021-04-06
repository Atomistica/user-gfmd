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
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


//extern __shared__ F_FLOAT sharedmem[]; // we will use this for reduce, when implemented



__global__ void Cuda_FixAddForce_DynamicCuda_PostForce_Kernel( // cu->cu_x, cu->cu_f,
                    int groupbit, // are these int? what do i use instead of const viscous doesnt need nlocal passed in
		    F_FLOAT reflxcenter_sigma, F_FLOAT reflycenter_sigma,
	            F_FLOAT radspersigma_x, F_FLOAT radspersigma_y,
		    F_FLOAT acontact_sigma, 
		    F_FLOAT border_dist,
		    F_FLOAT r_overEstar)//,
		    //F_FLOAT* foriginaladf)

{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;

  if (i < _nlocal) {
    if (_mask[i] & groupbit) {

      F_FLOAT arad_sq = pow(acontact_sigma, 2.0);
      F_FLOAT arad_border_sq = pow(acontact_sigma + border_dist, 2.0);

      F_FLOAT envelope, rad_sq;

      F_FLOAT xdist_to_center = _x[i + 0*_nmax] - reflxcenter_sigma;
      F_FLOAT ydist_to_center = _x[i + 1*_nmax] - reflycenter_sigma;

      // Look to closest periodic image (may be able to speed up on GPU)
      if (xdist_to_center > (_boxhi[0] - _boxlo[0])/2) xdist_to_center -= (_boxhi[0] - _boxlo[0]);
      if (xdist_to_center <-(_boxhi[0] - _boxlo[0])/2) xdist_to_center += (_boxhi[0] - _boxlo[0]);
      if (ydist_to_center > (_boxhi[1] - _boxlo[1])/2) ydist_to_center -= (_boxhi[1] - _boxlo[1]);
      if (ydist_to_center <-(_boxhi[1] - _boxlo[1])/2) ydist_to_center += (_boxhi[1] - _boxlo[1]);


      if (i < 4) printf("i%d adf xdisttocenter%f y%f arad_border_sq%f\n",i,xdist_to_center,ydist_to_center,arad_border_sq);//TAS
      //choice of corrugation
//#if 0      if (0) {
//        //potential is cropped sines to make cusps
//        cx = sin(fmod(xdist_to_center*radspersigma_x/2., M_PI) );
//        cy = sin(fmod(ydist_to_center*radspersigma_y/2., M_PI) );
//        sx = .5*cos(fmod(xdist_to_center*radspersigma_x/2., M_PI) );
//        sy = .5*cos(fmod(ydist_to_center*radspersigma_y/2., M_PI) );
//      } else {
//#endif
        // potential is cosines.  force is negative derivative (sines).
        // origin is the envelope center
        F_FLOAT cx= cos(xdist_to_center*radspersigma_x);
        F_FLOAT cy= cos(ydist_to_center*radspersigma_y);
        F_FLOAT sx= sin(xdist_to_center*radspersigma_x);
        F_FLOAT sy= sin(ydist_to_center*radspersigma_y);
//      }



      //if ((x[i][0] < 2) && (x[i][1] < 2)) 
      //  printf("x %f y %f cx cy sx sy %f %f %f %f\n", x[i][0], x[i][1], cx, cy, sx, sy);
      rad_sq = (pow(xdist_to_center, 2.0) + 
                pow(ydist_to_center, 2.0));

      if (i < 4) printf("i%d cx%f rad_sq%f\n",i,cx, rad_sq);//TAS



      // Choice of enevelope
//#if 0      if (0) {
//      // Envelope is a Hertzian profile
//        F_FLOAT p0Hertz=2.0 / M_PI * sqrt(aradsq) / r_overEstar; // fixed march 6 
//	F_FLOAT norm_to_latpot_coupling = 1; // assumed, but is probably nonlinear
//        if (rad_sq < aradsq) envelope = p0Hertz * sqrt(1.0 - rad_sq/aradsq);
//        else envelope = 0;
//      } else {
//#endif
      // Envelope is disk with height=1 and cosine roll-off over 1/2 nn spacing
        if (rad_sq < arad_sq) envelope = 1;
        else if (rad_sq < arad_border_sq) 
  	  envelope = 0.5 + 0.5*cos(M_PI*(sqrt(rad_sq)-acontact_sigma)/border_dist);
        else envelope = 0;
//      }
      // Envelope is diamond-tiled disk with height=1 so V=0 at edge... not implemented


      // foriginal[0] = "potential energy" for added force
      // foriginal[123] = force on atoms before extra force added
      //foriginaladf[0] += (double)(cx*envelope/radspersigma_x + cy*envelope/radspersigma_y); // added feb 12 2014
      //foriginaladf[1] += (double)(_f[i + 0 * _nmax]);
      //foriginaladf[2] += (double)(_f[i + 1 * _nmax]);
      //cuda_reduce(foriginal_loc,foriginal_total);

      if (i < 4) printf("i%d adfcontrib toxforce%f for total %f (cuz envelope %f)\n",i,sx*envelope, _f[i + 0 * _nmax]+ sx*envelope,envelope);//TAS
      if (i < 4) printf("i adfcontrib forces are at %p = %p \n",&(_f[i + 0 * _nmax]),&(_f[i]));//TAS

      _f[i + 0 * _nmax] += sx*envelope;
      _f[i + 1 * _nmax] += sy*envelope;

    }
  }
}

// #if 0 // This still has to be updated from the cpu version
// /* ---------------------------------------------------------------------- */
// 
// void FixAddForce_Dynamic::min_post_force(int vflag)
// {
//   //if (logfile) fprintf(logfile, "ADF calling post force from min_post_force\n");
//   post_force(vflag);
// }
// 
// /* ----------------------------------------------------------------------
//    potential energy of added force
// ------------------------------------------------------------------------- */
//  // There is a slight error in this calc due to the roll-off at the patch edge
// double FixAddForce_Dynamic::compute_scalar()
// {
//   // only sum across procs one time
// 
//   if (force_flag == 0) { // eflag is set to zero each time post_force called
//     MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
//     force_flag = 1;
//   }
//   return foriginal_all[0];
// }
// 
// /* ----------------------------------------------------------------------
//    potential energy of added force
// ------------------------------------------------------------------------- */
// 
// double FixAddForce_Dynamic::compute_vector(int n)
// {
//   // only sum across procs one time
// 
//   if (force_flag == 0) {
//     MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
//     force_flag = 1;
//   }
//   return foriginal_all[n+1];
// }
// #endif
