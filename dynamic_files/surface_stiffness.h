#ifndef __SURFACE_STIFFNESS_H
#define __SURFACE_STIFFNESS_H

#include "domain.h"
#include "error.h"
#include "force.h"
#include "linearalgebra.h"
#include "memory.h"

using namespace std;

class CrystalSurface;

namespace LAMMPS_NS {

class StiffnessKernel {
 public:
  StiffnessKernel(int, int *, char **, Domain *, Memory *, Error *);
  virtual ~StiffnessKernel();

  int get_dimension() {
    return dim_;
  };

  virtual const double *get_cell() const {
    return unit_cell;
  }

  virtual const double *get_invcell() const {
    return unit_cell;
  }

  int get_number_of_atoms() {
    return nu_;
  }

  int get_height() {
    return height_;
  }

  char *get_name() {
    return name_;
  };

  int get_number_of_parameters() {
    return npars_;
  }

  virtual bool get_mapping(int num_ref, double *ref_positions, int *indices,
                           int *types=NULL, double tol=1e-3) const {
    error_->all(FLERR,"This stiffness kernel cannot identify grid indices.");
    return false;
  }

  void set_invariant(bool invariant) {
    invariant_ = invariant;
  };

  /*!
   * Set model parameters
   */
  virtual void set_parameters(double *);

  /*!
   * Return 3x3 per atom dynamical matrices U, V, W, X, etc.
   */
  virtual void get_per_layer_dynamical_matrices(double, double,
						double_complex **,
						double *, double *);

  /*!
   * Called before get_dynamical_matrices
   */
  virtual void pre_compute() {
  };

  /*!
   * Called after get_dynamical_matrices
   */
  virtual void post_compute() { };

  /*!
   * Assemble per atom dynamical matrices into per cell matrices U and V
   */
  virtual void get_dynamical_matrices(double, double, double_complex *,
                                      double_complex *, double_complex *,
                                      double_complex dU=0.0);

  /*!
   * Renormalize interaction and compute surface stiffness
   */
  virtual void get_stiffness_matrix(double, double, double_complex *,
                                    double_complex dU=0.0);
  /*!
   * Renormalize interaction and compute surface stiffness for variable depth
   */
  virtual void get_stiffness_matrix(int, double, double, double_complex *,
				    double_complex dU=0.0);

  /*!
   * Renormalize interaction and compute surface stiffness - including grid
   * information. (Needed for all real-space methods.)
   */
  virtual void get_stiffness_matrix(int, int, double, double, double_complex *,
                                    double_complex dU=0.0);

  /*!
   * Get the n,0 Green's function matrix element
   */
  virtual void get_Gn0(int, double, double, double_complex *,
                       double_complex dU=0.0);

  /*!
   * Get linear force contribution at the gamma point (surface relaxation)
   */
  virtual void get_force_at_gamma_point(double *);

  /*!
   * Dump generic text info to file
   */
  virtual void dump_info(FILE *);

  /*!
   * Is this a manybody potential? This means that the energy of half of the
   * lattice planes will come from the potential and GFMD, respectively.
   */
  virtual bool is_manybody() { return false; }
  
 protected:
  /* Domain information, defines the Brillouin zone */
  Domain *domain_;

  /* Memory and error management */
  Memory *memory_;
  Error *error_;

  /* Dimension of the dynamical matrix at only nearest-layer interactions */
  int dim_;

  /* Number of degrees of freedom per unit cell */
  int ndof_;

  /* Number of layers for dynamical matrices with one atom per matrix */
  int nu_;

  /* Number of spring constants defined (1 for 1 layer fcc, 3 for 2 layers...) */
  int nk_;

  /* Number of layers in full sample */
  int height_;

  /* Number of (double) parameters for this model */
  int npars_;

  /* Should this kernel be translational invariant, i.e. observe q=0 sum
     rule */
  bool invariant_;

  /* Descriptor for this kernel */
  char name_[80];

  /* Stiffness kernel buffer */
  double_complex **D_;

  static double unit_cell[9];
};


StiffnessKernel *stiffness_kernel_factory(char *, int, int *, char **,
                                          Domain *, Force *, Memory *, Error *);


void iterate_Gnn(int, double_complex *, double_complex *, double_complex *,
		 double_complex *, double, int, Error *);
void iterate_Gn0(int, double_complex *, double_complex *, double_complex *,
		 double_complex *, double_complex *, int, Error *);
void direct_inversion_stiffness(int, int, int, double_complex *,
				double_complex *, double_complex *,
				double_complex *, Error *);
void displacement_transfer_matrix_stiffness(int, int, int, double_complex *,
					    double_complex *, double_complex *,
					    double_complex *, Error *);
void greens_function_transfer_matrix_stiffness(int, int, int, double_complex *,
					       double_complex *,
					       double_complex *,
					       double_complex *, Error *);
void transfer_matrix_Gn0(int, int, int, int,
                         double_complex *,
                         double_complex *,
                         double_complex *,
                         double_complex *,
                         double_complex *,
                         Error *);
void renormalization_group_stiffness(int, int, int, double_complex *,
				     double_complex *, double_complex *,
				     double_complex *, Error *);

void enforce_phi0_sum_rule(int, double_complex *, double);

}

#endif
