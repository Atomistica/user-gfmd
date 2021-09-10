#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "comm.h"
#include "pointers.h"
#include "time.h"
#ifndef NO_LAMMPS
#include "update.h"
#include "random_mars.h"
#endif

#include "gfmd_solver_dynamic.h"

/* ----------------------------------------------------------------------
 * Default dynamic implementation using the LAMMPS FFT wrapper
 * --------------------------------------------------------------------*/

#define round(x)  (int)(x < 0 ? (x - 0.5) : (x + 0.5))

GFMDSolverDynamic::GFMDSolverDynamic(LAMMPS *lmp, int narg, int *iarg,
                                     char **arg)
  : GFMDSolverFFT(lmp)
{
  strcpy(name, "dynamic");

  dump_ = NULL;

  if (*iarg > narg-5) {
    error->all(FLERR,"solver default/dynamic: Expected number of layers for q=0, number of layers for q=qmax, layer "
               "spacing, mass and dissipation arguments.");
  }

  char *endptr;

  if (*iarg >= narg)
    error->all(FLERR, "solver default/dynamic: Expected n(q=0) argument.");
  nmax_ = strtol(arg[*iarg], &endptr, 10);
  if (endptr == arg[*iarg]) {
    char errstr[1024];
    sprintf(errstr, "solver default/dynamic: Could not convert '%s' to "
            "integer.", arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;
  if (*iarg >= narg)
    error->all(FLERR, "solver default/dynamic: Expected n(q=qmax) argument.");
  nmin_ = strtol(arg[*iarg], &endptr, 10);
  // Min depth is 4 - with surface layer
  if (nmin_ <= 3) nmin_ = 3;
  if (endptr == arg[*iarg]) {
    char errstr[1024];
    sprintf(errstr, "solver default/dynamic: Could not convert '%s' to "
	    "integer.", arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;
  if (*iarg >= narg)
    error->all(FLERR, "solver default/dynamic: Expected delta argument.");
  delta_ = strtod(arg[*iarg], &endptr);
  if (endptr == arg[*iarg]) {
    char errstr[1024];
   sprintf(errstr, "solver default/dynamic: Could not convert '%s' to double.",
            arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;
  if (*iarg >= narg)
    error->all(FLERR, "solver default/dynamic: Expected mass argument.");
  mass_ = strtod(arg[*iarg], &endptr);
  if (endptr == arg[*iarg]) {
    char errstr[1024];
    sprintf(errstr, "solver default/dynamic: Could not convert '%s' to double.",
            arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;
  if (*iarg >= narg)
    error->all(FLERR, "solver default/dynamic: Expected gamma argument.");
  gamma_ = strtod(arg[*iarg], &endptr);
  if (endptr == arg[*iarg]) {
    char errstr[1024];
    sprintf(errstr, "solver default/dynamic: Could not convert '%s' to double.",
            arg[*iarg]);
    error->all(FLERR,errstr);
  }
  (*iarg)++;

  if (!comm->ghost_velocity){
    error->all(FLERR, "solver default/dynamic: Need to communicate ghost velocities.");
  }
  
  //  (*iarg)++;
  
  string bc_str;
  bc_flag_ = 0;
  if (*iarg < narg) {
    if (!strcmp(arg[*iarg], "bc")) {
      (*iarg)++;
      if (!strcmp(arg[*iarg], "rigid")){
	bc_flag_ = 1;
	printf("Using rigid boundary condition\n");
	bc_str = "using rigid boundary condition"; 
      }
      else if (!strcmp(arg[*iarg], "free")){
	bc_flag_ = 2;
	printf("Using free boundary condition\n");
	bc_str = "using free boundary condition";
      }
      else {
	bc_flag_ = 0;
	printf("Using phi boundary condition\n");
	bc_str = "using static GF boundary condition";
      }
      (*iarg)++;
    }
  }

  printf(" boundary condition string: %i\n",bc_flag_);
  (*iarg)++;


  lb_flag_ = 0;
  if (*iarg < narg) {
    if (!strcmp(arg[*iarg], "balance")) {
      lb_flag_ = 1;
      (*iarg)++;
    }
  }
  

  tflag_ = 0;
  /*
  if (*iarg < narg) {
    if (!strcmp(arg[*iarg], "langevin")) {
      tflag_ = 1;
      (*iarg)++;
      if (*iarg >= narg)
	error->all(FLERR, "solver default/dynamic: Expected temperature argument.");
      ttarget_ = strtod(arg[*iarg], &endptr);
      printf(" Langevin temp: %f\n",ttarget_);
      if (endptr == arg[*iarg]) {
	char errstr[1024];
	sprintf(errstr, "solver default/dynamic: Could not convert '%s' to double.",
		arg[*iarg]);
	error->all(FLERR,errstr);
      }
      (*iarg)++;
      if (*iarg >= narg)
        error->all(FLERR, "solver default/dynamic: Expected temperature seed argument.");
      tseed_ = strtol(arg[(*iarg)++], &endptr, 10);
      printf(" Langevin seed: %i\n",tseed_);
      if (tseed_ <= 0) error->all(FLERR,"Illegal fix langevin command");
      if (endptr == arg[*iarg]) {
        char errstr[1024];
        sprintf(errstr, "solver default/dynamic: Could not convert '%s' to int.",
                arg[*iarg]);
        error->all(FLERR,errstr);
      }
    }
    }*/

  if (*iarg < narg) {
    if (!strcmp(arg[*iarg], "dump")) {
      dump_ = fopen("gfmd_solver_dynamic.out", "w");
      (*iarg)++;
    }
  }    
    
  if (screen && comm->me == 0){
    fprintf(screen, "Using at least %i dynamic elastic layers spaced at %f for q>0.\n"
	    "Using %i dynamic elastic layers for q=0.\n"
	    "Mass of each atom within these layers is %f.\n"
	    "Dynamics is damped with a damping constant of %f*|q|.\n",
	    nmin_, delta_, nmax_, mass_, gamma_);
    if (tflag_)
      fprintf(screen,"Langevin dynamics are enabled with damping coefficient %f, target"
	      " temperature %f and initial seed %i.\n"
	      "Make sure to enable fix langevin for real-space dynamics.\n",
	      gamma_,ttarget_,tseed_);
  }
  if (logfile)
    fprintf(logfile, "Using at least %i dynamic elastic layers spaced at %f for q>0.\n"
	    "Using %i dynamic elastic layers for q=0.\n"
	    "Mass of each atom within these layers is %f.\n"
	    "Dynamics is damped with a damping constant of %f*|q|.\n",
	    nmin_, delta_, nmax_, mass_, gamma_);

  n_ = NULL;
  u_ = NULL;
  v_ = NULL;
  f_ = NULL;

  q_ = NULL;
  dyn_U0_ = NULL;
  dyn_U_ = NULL;
  dyn_V_ = NULL;
  dyn_G0_ = NULL;
  dyn_G_ = NULL;
  dyn_H_ = NULL;

  if (tflag_){
    random = NULL;
    random = new RanMars(lmp,tseed_ + comm->me);

    rn_ = NULL;
    
    dyn_sqrt_G0_ = NULL;
    dyn_sqrt_G_ = NULL;
  }

}


GFMDSolverDynamic::~GFMDSolverDynamic()
{
  if (dump_)
    fclose(dump_);

  if (tflag_){
    delete random;
    memory->destroy(dyn_sqrt_G0_);
    memory->destroy(dyn_sqrt_G_);

    if (n_) {
      for (int i = 0; i < nxy_loc; i++) {
	memory->destroy(rn_[i]);
      }

      delete rn_;
    }
  }
  
  if (n_) {
    for (int i = 0; i < nxy_loc; i++) {
      memory->destroy(u_[i]);
      memory->destroy(v_[i]);
      memory->destroy(f_[i]);
    }

    destroy_complex_buffer(u0_);
    destroy_complex_buffer(v0_);
    destroy_complex_buffer(f0_);

    delete n_;
    delete u_;
    delete v_;
    delete f_;
  }

  if (q_) {
    memory->destroy(q_);
    memory->destroy(dyn_U0_);
    memory->destroy(dyn_U_);
    memory->destroy(dyn_V_);
    memory->destroy(dyn_G0_);
    memory->destroy(dyn_G_);
    memory->destroy(dyn_H_);
  }
}


void GFMDSolverDynamic::set_grid_size(int in_nx, int in_ny, int in_ndof)
{
    GFMDSolverFFT::set_grid_size(in_nx, in_ny, in_ndof);

  n_ = new int[nxy_loc];
  u_ = new double_complex**[nxy_loc];
  v_ = new double_complex**[nxy_loc];
  f_ = new double_complex**[nxy_loc];
  
  u0_ = create_complex_buffer("GFMDSolverDynamic::u0");
  v0_ = create_complex_buffer("GFMDSolverDynamic::v0");
  f0_ = create_complex_buffer("GFMDSolverDynamic::f0");

  if (tflag_)
    rn_ = new double_complex**[nxy_loc];
    
  int nsum = 0;
  //nq_flag_ = 1;
  int idq = 0;

  int max_dim =  (nx>=ny) ? (nx) : (ny); 

  //printf("Proc ID: %i | x %i %i | y %i %i \n",comm->me,xlo_loc, xhi_loc,ylo_loc,yhi_loc);
  
  for (int i = xlo_loc; i <= xhi_loc; i++) {
    double qx = (i <= int((nx)/2)) ?
      (2.0*M_PI*(i)/nx) : (2.0*M_PI*(i-nx)/nx);
    for (int j = ylo_loc; j <= yhi_loc; j++) {
      double qy = (j <= int((ny)/2)) ?
        (2.0*M_PI*(j)/ny) : (2.0*M_PI*(j-ny)/ny);

      double qmag = sqrt(qx*qx+qy*qy);
      if( qmag < 0.000001) n_[idq] = nmax_;
      else {
	n_[idq] = round(nmax_/(qmag/(2.0*M_PI)*max_dim));
	if (n_[idq]<nmin_) n_[idq] = nmin_;
      }
      nsum += n_[idq];

      memory->create(u_[idq], n_[idq], ndof, "GFMDSolverDynamic::u");
      memory->create(v_[idq], n_[idq], ndof, "GFMDSolverDynamic::v");
      memory->create(f_[idq], n_[idq], ndof, "GFMDSolverDynamic::f");
      
      memset(u_[idq][0], 0, n_[idq]*ndof*sizeof(double_complex));
      memset(v_[idq][0], 0, n_[idq]*ndof*sizeof(double_complex));
      memset(f_[idq][0], 0, n_[idq]*ndof*sizeof(double_complex));


      if (tflag_){
	memory->create(rn_[idq], n_[idq]+1, ndof, "GFMDSolverDynamic::rn");
	memset(rn_[idq][0], 0, (n_[idq]+1)*ndof*sizeof(double_complex));
      }
      idq++;
    }
  }
  
  /*
  for (int i = 0; i < nxy_loc; i++) {
    if (i == 0 || nq_flag_ == 0){
      n_[i] = nmax_;
    }
    else {
      int ix = i/ny;
      int iy = i%ny;
      int ind_ix = (ix <= int((nx)/2)) ?
	(ix) : (nx-ix);
	int ind_iy = (iy<=int((ny)/2)) ?
	  (iy) : (ny-iy);
	
	n_[i] = round(nmax_/sqrt(ind_ix*ind_ix+ind_iy*ind_iy));
	if (n_[i]<nmin_) n_[i] = nmin_;
    }
    nsum += n_[i];
    memory->create(u_[i], n_[i], ndof, "GFMDSolverDynamic::u");
    memory->create(v_[i], n_[i], ndof, "GFMDSolverDynamic::v");
    memory->create(f_[i], n_[i], ndof, "GFMDSolverDynamic::f");
    
    memset(u_[i][0], 0, n_[i]*ndof*sizeof(double_complex));
    memset(v_[i][0], 0, n_[i]*ndof*sizeof(double_complex));
    memset(f_[i][0], 0, n_[i]*ndof*sizeof(double_complex));
  }
  
  double comp = 1.0*nsum/(nxy_loc*nmax_);
  if (screen && comm->me==0){ 
    fprintf(screen, "Total compression of layer count %f\n", comp); 
  }
  */
  printf("We ended up with nsum = %i \n",nsum);
  
}

void GFMDSolverDynamic::set_kernel(StiffnessKernel *kernel, bool normalize)
{
  memory->create(q_,      nxy_loc, "GFMDSolverDynamic::q");
  memory->create(dyn_U0_, nxy_loc, ndof_sq, "GFMDSolverDynamic::U0");
  memory->create(dyn_U_,  nxy_loc, ndof_sq, "GFMDSolverDynamic::U");
  memory->create(dyn_V_,  nxy_loc, ndof_sq, "GFMDSolverDynamic::V");
  memory->create(phi,     nxy_loc, ndof_sq, "GFMDSolverDynamic::phi");
  memory->create(dyn_G0_, nxy_loc, ndof_sq, "GFMDSolverDynamic::G0");
  memory->create(dyn_G_,  nxy_loc, ndof_sq, "GFMDSolverDynamic::G");
  memory->create(dyn_H_,  nxy_loc, ndof_sq, "GFMDSolverDynamic::H");

  if (tflag_){
    memory->create(dyn_sqrt_G0_, nxy_loc, ndof_sq, "GFMDSolverDynamic::sqrtG0");
    memory->create(dyn_sqrt_G_,  nxy_loc, ndof_sq, "GFMDSolverDynamic::sqrtG");
  }
  
  double inv_nxny = 1./double(nx*ny);
  
  int idq = 0;
  int n_static;
  double cx, cy, cx2, cy2, cxy2, sxy2;
  double_complex shift_factor;
  double_complex Gvar, G0var;

  double g0_r, g_r, g0_i, g_i, g0_mag, g_mag;
  
  for (int i = xlo_loc; i <= xhi_loc; i++) {
    double qx = (i <= int((nx)/2)) ?
      (2.0*M_PI*(i)/nx) : (2.0*M_PI*(i-nx)/nx);
    for (int j = ylo_loc; j <= yhi_loc; j++) {
      double qy = (j <= int((ny)/2)) ? 
        (2.0*M_PI*(j)/ny) : (2.0*M_PI*(j-ny)/ny);

      memset(dyn_G0_[idq], 0, ndof_sq*sizeof(double_complex));
      memset(dyn_G_[idq],  0, ndof_sq*sizeof(double_complex));
      memset(dyn_H_[idq],  0, ndof_sq*sizeof(double_complex));

      if (tflag_){
	memset(dyn_sqrt_G0_[idq], 0, ndof_sq*sizeof(double_complex));
	memset(dyn_sqrt_G_[idq],  0, ndof_sq*sizeof(double_complex));
      }
      
      q_[idq] = sqrt( qx*qx + qy*qy );
      kernel->get_dynamical_matrices(qx, qy, dyn_U0_[idq], dyn_U_[idq],
                                     dyn_V_[idq]);

      // Define damping matrices
      cx = cos(qx);
      cy = cos(qy);
      cx2 = cos(qx/2.0);
      cy2 = cos(qy/2.0);

      cxy2 = cos( (qx+qy)/2.0 );
      sxy2 = sin( (qx+qy)/2.0 );

      shift_factor = cexp(COMPLEX_NUMBER(0.0, 0.5*(qx+qy)));

      //printf("cx %f cy %f cx2 %f cy2 %f cxy2 %f sxy2 %f\n",cx,cy,cx2,cy2,cxy2,sxy2);
      //printf("%f %f %f %f %f %f\n",qx, qy, creal(factor), cimag(factor), creal(conj(factor)), cimag(conj(factor)));
            
      MEL(ndof, dyn_G0_[idq], 0, 0) = gamma_ * (8 - 2.0*(cx + cy));
      MEL(ndof, dyn_G0_[idq], 1, 1) = gamma_ * (8 - 2.0*(cx + cy));
      MEL(ndof, dyn_G0_[idq], 2, 2) = gamma_ * (8 - 2.0*(cx + cy));

      MEL(ndof, dyn_G_[idq], 0, 0) = gamma_ * (12 - 2.0*(cx + cy));
      MEL(ndof, dyn_G_[idq], 1, 1) = gamma_ * (12 - 2.0*(cx + cy));
      MEL(ndof, dyn_G_[idq], 2, 2) = gamma_ * (12 - 2.0*(cx + cy));

      MEL(ndof, dyn_H_[idq], 0, 0) = -4 * gamma_ * cx2 * cy2 * shift_factor;
      MEL(ndof, dyn_H_[idq], 1, 1) = -4 * gamma_ * cx2 * cy2 * shift_factor;
      MEL(ndof, dyn_H_[idq], 2, 2) = -4 * gamma_ * cx2 * cy2 * shift_factor;

      if (tflag_){
	/*
	MEL(ndof, dyn_sqrt_G0_[idq], 0, 0) = sqrt(gamma_ * (8 - 2.0*(cx + cy)));
	MEL(ndof, dyn_sqrt_G0_[idq], 1, 1) = sqrt(gamma_ * (8 - 2.0*(cx + cy)));
	MEL(ndof, dyn_sqrt_G0_[idq], 2, 2) = sqrt(gamma_ * (8 - 2.0*(cx + cy)));
	
	MEL(ndof, dyn_sqrt_G_[idq], 0, 0) = sqrt(gamma_ * (12 - 2.0*(cx + cy)));
	MEL(ndof, dyn_sqrt_G_[idq], 1, 1) = sqrt(gamma_ * (12 - 2.0*(cx + cy)));
	MEL(ndof, dyn_sqrt_G_[idq], 2, 2) = sqrt(gamma_ * (12 - 2.0*(cx + cy)));
	*/

	//Gvar =  12 - 2.0*(cx + cy) - 8 * cx2 * cy2 * shift_factor * cxy2 ;
	//G0var =  8 - 2.0*(cx + cy) - 4 * cx2 * cy2 * shift_factor * shift_factor ;

	//g_r = 12 - 2 * (cx + cy) - 8 * cx2 * cy2 * cxy2 * cxy2;
	//g_i = 8 * cx2 * cy2 * sxy2 * cxy2;

	//g_mag = sqrt(g_r * g_r + g_i * g_i);
	
	g0_r = 8 - 2.0 * (cx + cy) - 4.0 * cx2 * cy2 * cxy2;
        g0_i = -4.0 * cx2 * cy2 * sxy2;

        g0_mag = sqrt(g0_r * g0_r + g0_i * g0_i);


	MEL(ndof, dyn_sqrt_G0_[idq], 0, 0) = COMPLEX_NUMBER( sqrt( 0.5 * gamma_ * (g0_mag + g0_r) ), sqrt( 0.5 * gamma_ * (g0_mag - g0_r) ) );
	MEL(ndof, dyn_sqrt_G0_[idq], 1, 1) = COMPLEX_NUMBER( sqrt( 0.5 * gamma_ * (g0_mag + g0_r) ), sqrt( 0.5 * gamma_ * (g0_mag - g0_r) ) );
	MEL(ndof, dyn_sqrt_G0_[idq], 2, 2) = COMPLEX_NUMBER( sqrt( 0.5 * gamma_ * (g0_mag + g0_r) ), sqrt( 0.5 * gamma_ * (g0_mag - g0_r) ) );
		
	MEL(ndof, dyn_sqrt_G_[idq], 0, 0) = COMPLEX_NUMBER( sqrt( gamma_ * ( 12 - 2.0*(cx + cy) - 8.0 * cx2 * cy2 * cxy2 )), 0);
	MEL(ndof, dyn_sqrt_G_[idq], 1, 1) = COMPLEX_NUMBER( sqrt( gamma_ * ( 12 - 2.0*(cx + cy) - 8.0 * cx2 * cy2 * cxy2 )), 0);
        MEL(ndof, dyn_sqrt_G_[idq], 2, 2) = COMPLEX_NUMBER( sqrt( gamma_ * ( 12 - 2.0*(cx + cy) - 8.0 * cx2 * cy2 * cxy2 )), 0);
	
	/*
	printf("idq %i qx %f qy %f G %f %f %f %f %f %f\n",idq,qx,qy,creal(dyn_sqrt_G0_[idq][0]),
	       cimag(dyn_sqrt_G0_[idq][0]), creal(dyn_sqrt_G0_[idq][4]), cimag(dyn_sqrt_G0_[idq][4]),
	       creal(dyn_sqrt_G0_[idq][8]), cimag(dyn_sqrt_G0_[idq][8]));
	*/
	
      }
      
      /*
      printf("idq %i qx %f qy %f G %f %f %f %f %f %f %f %f %f\n",idq,qx,qy,creal(dyn_G0_[idq][0]),
             creal(dyn_G0_[idq][1]), creal(dyn_G0_[idq][2]), creal(dyn_G0_[idq][3]),
	     creal(dyn_G0_[idq][4]), creal(dyn_G0_[idq][5]),
             creal(dyn_G0_[idq][6]), creal(dyn_G0_[idq][7]), creal(dyn_G0_[idq][8]));      
      printf("idq %i qx %f qy %f G %f %f %f %f %f %f %f %f %f\n\n",idq,qx,qy,sqrt(creal(dyn_G0_[idq][0])),
	     sqrt(creal(dyn_G0_[idq][1])),sqrt(creal(dyn_G0_[idq][2])),sqrt(creal(dyn_G0_[idq][3])),
	     sqrt(creal(dyn_G0_[idq][4])),sqrt(creal(dyn_G0_[idq][5])),sqrt(creal(dyn_G0_[idq][6])),
	     sqrt(creal(dyn_G0_[idq][7])),sqrt(creal(dyn_G0_[idq][8])));
      */
      // Dynamic + static depth
      n_static = nmax_ - n_[idq] ;
      kernel->get_stiffness_matrix( n_static, qx, qy, phi[idq]);

      // Last layer has full neighbors
      for (int k=0; k<ndof_sq; k++){
	phi[idq][k] += dyn_U_[idq][k] - dyn_U0_[idq][k];
      }

      idq++;
    }
  }
}


double GFMDSolverDynamic::post_force(void *input_buffer_ptr,
				     void *output_buffer_ptr,
				     char *dump_prefix)
{
  double ft_time = 0;
  double it_time = 0;

  clock_t start;

  start = clock();
  double **input_buffer = static_cast<double**>(input_buffer_ptr);
  double **output_buffer = static_cast<double**>(output_buffer_ptr);
  
  fft_forward(input_buffer, u0_);
  fft_forward(output_buffer, v0_);

  ft_time += 1.0*(clock()-start)/(double) CLOCKS_PER_SEC;
  
  if (dump_prefix)
    dump(dump_prefix, u0_);

  start = clock();
  verlet_step1();
  double epot = energy_and_forces();
  if (tflag_) langevin_thermostat();
  double ekin = verlet_step2();

  it_time = 1.0*(clock()-start)/(double) CLOCKS_PER_SEC;

  if (dump_)
    fprintf(dump_, "%e %e  %e %e\n", epot0_, epot1_, ekin0_, ekin1_);

  start = clock();
  fft_reverse(f0_, output_buffer);
  ft_time += 1.0*(clock()-start)/(double) CLOCKS_PER_SEC;
  
  //printf("xx %i %e %e\n",comm->me, ft_time, it_time );
  return epot + ekin;
}


double GFMDSolverDynamic::energy_and_forces()
{

  double epot0 = 0.0, epot1 = 0.0;
  memset(f0_[0], 0, nxy_loc*ndof*sizeof(double_complex));

  double inv_nxny = 1./double(nx*ny);

  for (int idq = 0; idq < nxy_loc; idq++) {

    int ni = n_[idq];
    double q2 = q_[idq]*q_[idq];

    double epot = 0.0;
    double_complex *ui0 = u0_[idq];
    double_complex *vi0 = v0_[idq];
    double_complex *fi0 = f0_[idq];
    double_complex **ui = u_[idq];
    double_complex **vi = v_[idq];
    double_complex **fi = f_[idq];

    double_complex *U0 = dyn_U0_[idq];
    double_complex *U  = dyn_U_[idq];
    double_complex *V  = dyn_V_[idq];

    double_complex *G0 = dyn_G0_[idq];
    double_complex *G  = dyn_G_[idq];
    double_complex *H  = dyn_H_[idq];

    double_complex *phi_q = phi[idq];
    
    for (int j = 0; j < ndof ; j++){
      ui0[j] *= inv_nxny;
      vi0[j] *= inv_nxny;
    }

    memset(fi[0], 0, ni*ndof*sizeof(double_complex));

    // Surface layer
    mat_mul_sub_vec(ndof, U0, ui0,   fi0);
    mat_mul_sub_vec(ndof, V,  ui[0], fi0);
    epot -= creal(conj_vec_dot_vec(ndof, ui0, fi0));

    mat_mul_sub_vec(ndof, G0, vi0,   fi0);
    mat_mul_sub_vec(ndof, H,  vi[0], fi0);

    // First GF layer
    conj_mat_mul_sub_vec(ndof, V, ui0,   fi[0]);
    mat_mul_sub_vec     (ndof, U, ui[0], fi[0]);
    mat_mul_sub_vec     (ndof, V, ui[1], fi[0]);
    epot -= creal(conj_vec_dot_vec(ndof, ui[0], fi[0]));

    conj_mat_mul_sub_vec(ndof, H, vi0,   fi[0]);
    mat_mul_sub_vec     (ndof, G, vi[0], fi[0]);
    mat_mul_sub_vec     (ndof, H, vi[1], fi[0]);

    // Intermediate layers
      for (int idn = 1; idn < ni-1; idn++) {
	conj_mat_mul_sub_vec(ndof, V, ui[idn-1], fi[idn]);
	mat_mul_sub_vec     (ndof, U, ui[idn],   fi[idn]);
	mat_mul_sub_vec     (ndof, V, ui[idn+1], fi[idn]);
	epot -= creal(conj_vec_dot_vec(ndof, ui[idn], fi[idn]));
	
	conj_mat_mul_sub_vec(ndof, H, vi[idn-1], fi[idn]);
	mat_mul_sub_vec     (ndof, G, vi[idn],   fi[idn]);
	mat_mul_sub_vec     (ndof, H, vi[idn+1], fi[idn]);
      }
    
    // Bottom layer
    conj_mat_mul_sub_vec(ndof, V,     ui[ni-2], fi[ni-1]);

    //bc_flag_ = 0;
    if (bc_flag_ == 0){
      mat_mul_sub_vec   (ndof, phi_q, ui[ni-1], fi[ni-1]);
      epot -= creal(conj_vec_dot_vec(ndof, ui[ni-1], fi[ni-1]));

      conj_mat_mul_sub_vec(ndof, H, vi[ni-2], fi[ni-1]);
      mat_mul_sub_vec     (ndof, G, vi[ni-1], fi[ni-1]);
    }
    
    else if (bc_flag_ == 1){
      mat_mul_sub_vec   (ndof, U,     ui[ni-1], fi[ni-1]);
      epot -= creal(conj_vec_dot_vec(ndof, ui[ni-1], fi[ni-1]));

      conj_mat_mul_sub_vec(ndof, H, vi[ni-2], fi[ni-1]);
      mat_mul_sub_vec     (ndof, G, vi[ni-1], fi[ni-1]);
    }

    else if (bc_flag_ == 2){
      mat_mul_sub_vec   (ndof, U0,     ui[ni-1], fi[ni-1]);
      epot -= creal(conj_vec_dot_vec(ndof, ui[ni-1], fi[ni-1]));

      conj_mat_mul_sub_vec(ndof, H, vi[ni-2], fi[ni-1]);
      mat_mul_sub_vec     (ndof, G0, vi[ni-1], fi[ni-1]);
    }

    if (idq == 0) {
      epot0 += epot;
    }
    else
      epot1 += epot;
   }
  
  epot0_ = 0.5*epot0;
  epot1_ = 0.5*epot1;
  return (epot0+epot1);
}

void GFMDSolverDynamic::langevin_thermostat()
{
  /*  
  // Correct for adding variances of 2 RN
  double inv_sqrt_2 = 1.0/sqrt(2.);
  lgprefactor = inv_sqrt_2*sqrt(mass_*24.0*force->boltz/update->dt/force->mvv2e*ttarget_) / force->ftm2v;
  //lgprefactor = inv_sqrt_2*sqrt(mass_*2.0*force->boltz/update->dt/force->mvv2e*ttarget_) / force->ftm2v;
  
  // Zero out random numbers
  for (int idq = 1; idq < nxy_loc; idq++) {
    memset(rn_[idq][0], 0, (n_[idq]+1)*ndof*sizeof(double_complex));
  }

  // Enforce symmetry of FT
  int nxy = nx*ny;
  for (int idq = 1; idq < nxy_loc; idq++) {
    if ((idq%ny)==0){
      idqm = nxy - idq;
    }
    else {
      idqm = nxy + nx - idq;
      idqm -= nxy*(idqm/nxy);
    }

    double_complex **rn1 = rn_[idq];
    double_complex **rn2 = rn_[idqm];
    int ni = n_[idq];
    
    // ni+1 layers
    double tmp;
    for (int idn = 0; idn <= ni; idn++) {
      for (int idof = 0; idof < ndof; idof++) {	
	tmp = lgprefactor*(random->uniform()-0.5);
	//tmp = lgprefactor*random->gaussian();
	
	rn1[idn][idof] += COMPLEX_NUMBER(tmp,0);
	rn2[idn][idof] += COMPLEX_NUMBER(tmp,0);

      }
    }
  }
  */

  // In general, might need to change units for gamma
  lgprefactor = sqrt(mass_*24.0*force->boltz/update->dt/force->mvv2e*ttarget_) / force->ftm2v;

  /*
  // Zero out random numbers
  for (int idq = 1; idq < nxy_loc; idq++) {
    memset(rn_[idq][0], 0, (n_[idq]+1)*ndof*sizeof(double_complex));
  }
  
  for (int idq = 1; idq < nxy_loc; idq++) {
    double_complex **rn1 = rn_[idq];
    int ni = n_[idq];
    double tmp;
    for (int idn = 0; idn <= ni; idn++) {
      for (int idof = 0; idof < ndof; idof++) {
	tmp = lgprefactor*(random->uniform()-0.5);
	rn1[idn][idof] = COMPLEX_NUMBER(tmp,0);
      }
    }
  }
  */
  
  // Skip q = 0
  for (int idq = 1; idq < nxy_loc; idq++) {
    
    double_complex *fi0 = f0_[idq];
    double_complex **fi = f_[idq];
    double_complex **rni = rn_[idq];
        
    int ni = n_[idq];

    double_complex *sqrtG0 = dyn_sqrt_G0_[idq];
    double_complex *sqrtG  = dyn_sqrt_G_[idq];

    //printmat(3,sqrtG0);
    //printf("\n");
    //printvec(3,fi0);
    //printf("\n");
    // Top layer
    //printvec(3,rni[0]);
    //printf("\n");

    //mat_mul_add_vec(ndof, sqrtG0, rni[0], fi0);

    //for (int idof = 0; idof < ndof; idof++) {
    //  fi0[idof] -= sqrtG0[4*idof]*lgprefactor*(random->uniform()-0.5);
    //}
    
    //printvec(3,fi0);
    //printf("\n\n");
    
    // Intermediate layers
    for (int idn = 0; idn < ni-1; idn++) {
      //printf("\n");
      //printvec(3,fi[idn]);
      //mat_mul_add_vec(ndof, sqrtG, rni[idn+1], fi[idn]);
      for (int idof = 0; idof < ndof; idof++) {
	fi[idn][idof] -= sqrtG[4*idof]*lgprefactor*(random->uniform()-0.5);
      }
      //printvec(3,fi[idn]);
      //printf("\n");
    }

    // Bottom layer
    bc_flag_ == 0;
    if (bc_flag_ == 2){
      //mat_mul_add_vec(ndof, sqrtG0, rni[ni], fi[ni-1]);
      for (int idof = 0; idof < ndof; idof++) {
        fi[ni-1][idof] -= sqrtG0[4*idof]*lgprefactor*(random->uniform()-0.5);
      }
    }
    else {
      //mat_mul_add_vec(ndof, sqrtG, rni[ni], fi[ni-1]);
      for (int idof = 0; idof < ndof; idof++) {
        fi[ni-1][idof] -= sqrtG[4*idof]*lgprefactor*(random->uniform()-0.5);
      }
    }

  }
}

void GFMDSolverDynamic::verlet_step1()
{
#ifndef NO_LAMMPS
  double dt = update->dt;
  dt *= force->ftm2v;
#else
  double dt = 0.1; 
#endif
  double d2t = update->dt*dt;
  
  for (int idq = 0; idq < nxy_loc; idq++) {
    double_complex **ui = u_[idq];
    double_complex **vi = v_[idq];
    double_complex **fi = f_[idq];

    for (int idn = 0; idn < n_[idq]; idn++) {
      for (int idof = 0; idof < ndof; idof++) {
	ui[idn][idof] += vi[idn][idof] * dt + 0.5 * fi[idn][idof] / mass_ * d2t;
	vi[idn][idof] +=                      0.5 * fi[idn][idof] / mass_ * dt;

      }
    }
  }
}


double GFMDSolverDynamic::verlet_step2()
{

#ifndef NO_LAMMPS
  double dt = update->dt;
  dt *= force->ftm2v;
#else
  double dt = 0.1;
#endif
  double d2t = update->dt*dt;
  double ekin0 = 0.0, ekin1 = 0.0;

  // Contributions to KE: 3(nxy-1)nz for q > 0
  for (int idq = 0; idq < nxy_loc; idq++) {
    double ekin = 0.0;
    double_complex **vi = v_[idq];
    double_complex **fi = f_[idq];

    for (int idn = 0; idn < n_[idq]; idn++) {
      for (int idof = 0; idof < ndof; idof++) {
	vi[idn][idof] += 0.5 * fi[idn][idof] / mass_ * dt;
	ekin += 0.5*mass_*creal(conj(vi[idn][idof])*vi[idn][idof]);
      }
          
    }
    if (idq == 0){
      ekin0 += ekin;
    }
    else
      ekin1 += ekin;
  }

  ekin0_ = ekin0;
  ekin1_ = ekin1;

  return ekin0+ekin1;
}


double GFMDSolverDynamic::memory_usage()
{
  double bytes = 0.0;

  // u, v, f
  for (int i = 0; i < nxy_loc; i++) {
    bytes += 3*n_[i]*ndof*sizeof(double_complex);
  }
  // dyn_U0, dyn_U, dyn_V, dyn_G0_, dyn_G_, dyn_H_
  bytes += 6*nxy_loc*ndof_sq*sizeof(double_complex);

  return bytes + GFMDSolverFFT::memory_usage();
}
