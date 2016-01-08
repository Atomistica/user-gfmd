#ifndef CONTACT_H
#define CONTACT_H

#include <fftw3.h>

#define NDOF 1
#define NDOF_SQ (NDOF*NDOF)

#include "pointers.h"
#include "linearalgebra.h"

namespace LAMMPS_NS {

class Contact: protected Pointers {
 public:
  Contact(const char *);
  ~Contact();

  double energy_and_forces(bool[NDOF], double);
  void compute(bool[NDOF], double);
  int min_tr(bool[NDOF], double);
  void loop();

  void get_f(double &, double &, double &);

 private:
  class FixContactMap *map;        // contact map
  class StiffnessKernel *kernel;   // stiffness kernel
#ifdef HAVE_SDL
  class LiveView *live_view;
#endif

  int nx_, ny_, nxy_, nq_;
  int max_iter, nrep, natt, ntol, nrepprevious, nattprevious;
  int nsteps_, outevery, outcounter;
  int mode_;

  char dump_fn[1024];
  char *restart_fn;

  double epot0_el_, epot_el_, epot_wall_, epot_ext_;
  double curv_wall[NDOF_SQ], ftol, rho_tol, dmax, atol;
  double z_, z0_, z1_, frep, fatt, delta, h0_, avgh_;

  double **u_r_, **f_r_, *fftr_;
  double_complex **u_q_, **f_q_, *fftq_;
  double_complex **phi_;

  bool fix_center_[NDOF];

  fftw_plan *fft_forward_u, *fft_reverse_u;
  fftw_plan *fft_forward_f, *fft_reverse_f;

  void ur_to_uq();
  void uq_to_ur();
  void fr_to_fq();
  void fq_to_fr();
  void xr_to_uq();
  void uq_to_xr();

  double get_fmaxq();
  double get_fmaxr();

  double get_fnorm();

  void wall_interaction();

  void print_status(char, int, double, double, double, double, double, int,
		    int, double, double);

};

}

#endif