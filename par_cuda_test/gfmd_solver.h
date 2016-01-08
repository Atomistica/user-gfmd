#ifndef GFMD_SOLVER_H
#define GFMD_SOLVER_H

#ifdef GFMD_CUFFT
#include <cufft.h>
#endif

#include "pointers.h"
#include "fft3d_wrap.h"

#include "linearalgebra.h"

namespace LAMMPS_NS {

class GFMDSolver : protected Pointers {
 public:
  GFMDSolver(LAMMPS *, size_t, size_t, size_t);
  virtual ~GFMDSolver();

  virtual void set_operator(double_complex **) = 0;

  // the solver is called asynchronously
  virtual void detach(double **, double **) { };
  virtual double join(double **, double **) = 0;

  double **create_double_buffer(char *);
  void destroy_double_buffer(double **);

  double_complex **create_complex_buffer(char *);
  void destroy_complex_buffer(double_complex **);

  double_complex **create_complex_operator_buffer(char *);
  void destroy_complex_operator_buffer(double_complex **);

  size_t get_xlo_loc() {
    return xlo_loc;
  }
  size_t get_xhi_loc() {
    return xhi_loc;
  }
  size_t get_ylo_loc() {
    return ylo_loc;
  }
  size_t get_yhi_loc() {
    return yhi_loc;
  }
  size_t get_nxy_loc() {
    return nxy_loc;
  }

  char *get_name() {
    return name;
  };

  virtual double memory_usage();

 protected:
  /* Descriptor for this grid class */
  char name[80];

  /* Parameters */
  size_t nx, ny, ndof, ndof_sq;

  /* Derived FFT stuff */

  /* my process rank, total number of procs */
  int me, nprocs;

  /* job assignment info */
  size_t xlo_loc, xhi_loc, nx_loc;
  size_t ylo_loc, yhi_loc, ny_loc;
  size_t nxy_loc;
};



class GFMDSolverDefault : public GFMDSolver {
 public:
  GFMDSolverDefault(LAMMPS *, size_t, size_t, size_t);
  ~GFMDSolverDefault();

  void set_operator(double_complex **);
  double join(double **, double **);

  double memory_usage();

 protected:
  FFT3d *fft;
  double *fft_data;
  double_complex **q_buffer;
  double_complex **q_operator;
};



#ifdef GFMD_CUFFT
class CUDAThread : public GFMDSolver {
 public:
  double **input_buffer, **output_buffer;
  double epot;

  CUDAThread(LAMMPS *, size_t, size_t, size_t, double_complex **);
  ~CUDAThread();

  void init();
  void run();
  double apply_operator();
  void del();

  void lock() { pthread_mutex_lock(&mutex); }
  void unlock() { pthread_mutex_unlock(&mutex); }
  void wait() { pthread_cond_wait(&cond, &mutex); }
  void signal() { pthread_cond_signal(&cond); }
  void stop() { thread_is_running = false; signal(); unlock(); }

  void set_operator(double_complex **);
  double join(double **, double **) { }

 protected:
  bool thread_is_running;
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  int cuda_device_id;

  size_t nbuffer_half_loc;

  cufftHandle fft_plan_forward, fft_plan_inverse;

  double_complex **q_operator;

  cufftDoubleReal *real_buffer;
  cufftDoubleComplex *complex_buffer;

  cufftDoubleComplex *dev_q_operator;
  size_t dev_q_operator_pitch;

  cufftDoubleReal *dev_u, *dev_f;
  size_t dev_u_pitch, dev_f_pitch;

  cufftDoubleComplex *dev_q_buffer1, *dev_q_buffer2;
  size_t dev_q_buffer_pitch;

  cufftDoubleReal *dev_epot_buf; 
};
#endif

class GFMDSolverCUDA : public GFMDSolver {
 public:
  GFMDSolverCUDA(LAMMPS *, size_t, size_t, size_t);
  ~GFMDSolverCUDA();

  void set_operator(double_complex **);
  void detach(double **, double **);
  double join(double **, double **);

#ifdef GFMD_CUFFT
 protected:
  pthread_t pthread;
  CUDAThread *thread;

  static void *cuda_thread_dispatcher(void *);
#endif
};



GFMDSolver *gfmd_solver_factory(char *, LAMMPS *, size_t, size_t, size_t);

}

#endif
