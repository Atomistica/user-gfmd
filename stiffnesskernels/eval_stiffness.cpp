/*
 * Standalone code, compute stiffness matrix within the first Brillouin
 * zone.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "error.h"
#include "memory.h"

#include "gfmd_stiffness.h"


struct spec_pt_t {
  char sym;
  double x;
  double y;
};


spec_pt_t spec_pts[] = {
  { 'G',  0.0,  0.0 },
  { 'M',  0.5,  0.0 },
  { 'R',  0.5,  0.5 },
  { 'm', -0.5,  0.0 },
  { 'r', -0.5,  0.5 },
  { 0,    0.0,  0.0 }
};


void find_spec_pt_coord(char sym, double *x, double *y)
{
  int i;

  *x = -10.0;
  *y = -10.0;

  i = 0;
  while (spec_pts[i].sym != 0) {
    if (spec_pts[i].sym == sym) {
      *x = spec_pts[i].x;
      *y = spec_pts[i].y;
    }
    i++;
  }

  if (*x < -1.0) {
    printf("\nWarning: Could not find special point %c.\n", sym);
    *x = 0.0;
    *y = 0.0;
  }
}


void band_lines(int nspec, double *specx, double *specy, double *specd,
		int npts, double *ptx, double *pty, double *ptd)
{
  int i, j, k, ncpts;
  double l, cl, kl, dl, x0, y0, x1, y1, d;

  l = 0.0;
  x0 = specx[0];
  y0 = specy[0];
  for (i = 1; i < nspec; i++) {
    x1 = specx[i];
    y1 = specy[i];
    l += sqrt(pow(x1-x0, 2) + pow(y1-y0, 2));
    x0 = x1;
    y0 = y1;
  }

  d  = (npts-1)/l;
  dl = 1.0/d;

  k  = 0;
  kl = 0.0;
  x0 = specx[0];
  y0 = specy[0];
  specd[0] = 0.0;
  for (i = 1; i < nspec; i++) {
    x1 = specx[i];
    y1 = specy[i];
    cl = sqrt(pow(x1-x0, 2) + pow(y1-y0, 2));

    if (i == nspec-1) {
      ncpts = npts - k - 1;
    }
    else {
      ncpts = (int) (d*cl);
    }

    for (j = 0; j < ncpts; j++) {
      ptx[k] = x0 + (x1-x0)*j/ncpts;
      pty[k] = y0 + (y1-y0)*j/ncpts;
      ptd[k] = kl;
      k++;
      kl += dl;
    }

    x0 = x1;
    y0 = y1;

    specd[i] = kl;
  }

  ptx[npts-1] = specx[nspec-1];
  pty[npts-1] = specy[nspec-1];
  ptd[npts-1] = l;
}


#define MAX_DYN_DIM_SQ 1024
#define SYM_LEN 32

int main(int narg, char *arg[])
{
  int iarg, nspec, npts, i;
  double *specx, *specy, *specd;
  double *ptx, *pty, *ptd;
  double_complex mat[MAX_DYN_DIM_SQ];
  char *specsym, *stifffn, *gffn;

  StiffnessKernel *kernel;

  Domain domain;
  Error error;
  Memory memory(&error);

  FILE *f;

  /* --- */

  domain.xprd = 1.0;
  domain.yprd = 1.0;
  domain.zprd = 1.0;

  if (narg < 4) {
    printf("Syntax: eval_stiffness <filename with special points> "
	   "<stiffness filename> <greens functions filename> "
	   "<number of band points> <kernel> <kernel parameters>\n");
    exit(999);
  }

  /*
  nspec = atoi(arg[1]);
  specsym = (char *) malloc(nspec*sizeof(char));
  specx = (double *) malloc(nspec*sizeof(double));
  specy = (double *) malloc(nspec*sizeof(double));
  specd = (double *) malloc(nspec*sizeof(double));
  */

  iarg = 1;

  f = fopen(arg[iarg], "r");
  if (!f) {
    printf("Error opening %s.\n", arg[iarg]);
    exit(999);
  }
  iarg++;
  fscanf(f, "%i", &nspec);
  specsym = (char *) malloc(nspec*SYM_LEN*sizeof(char));
  specx = (double *) malloc(nspec*sizeof(double));
  specy = (double *) malloc(nspec*sizeof(double));
  specd = (double *) malloc(nspec*sizeof(double));
  for (i = 0; i < nspec; i++) {
    fscanf(f, "%s%lf%lf", &specsym[i*SYM_LEN], &specx[i], &specy[i]);
  }
  fclose(f);

  stifffn = arg[iarg];
  iarg++;
  gffn = arg[iarg];
  iarg++;

  npts = atoi(arg[iarg]);
  iarg++;

  if (npts <= nspec) {
    printf("Number of band points must be larger that number of special "
	   "points.\n");
    exit(999);
  }

  iarg++;
  kernel = stiffness_kernel_factory(arg[iarg-1], narg, &iarg, arg,
				    &domain, &memory, &error);
  if (!kernel) {
    printf("Could not find stiffness kernel %s.\n", arg[iarg-1]);
    exit(999);
  }

  ptx = (double *) malloc(npts*sizeof(double));
  pty = (double *) malloc(npts*sizeof(double));
  ptd = (double *) malloc(npts*sizeof(double));

  band_lines(nspec, specx, specy, specd, npts, ptx, pty, ptd);

  /*
   * stiffness.out
   */
  f = fopen(stifffn, "w");
  if (!f) {
    printf("Error opening %s.\n", stifffn);
    exit(999);
  }
  fprintf(f, "# dist  qx qy");
  int dim = kernel->get_dimension();
  i = 4;
  for (int k = 0; k < dim; k++)
    for (int l = 0; l < dim-k; l++) {
      fprintf(f, "  %i:mat%i%i", i, l+1, l+k+1);
      i += 2;
    }
  fprintf(f, "\n");
  for (i = 0; i < npts; i++) {
    kernel->get_stiffness_matrix(2*M_PI*ptx[i], 2*M_PI*pty[i], mat);
    fprintf(f, "%e  %e %e", ptd[i], ptx[i], pty[i]);
    for (int k = 0; k < dim; k++)
      for (int l = 0; l < dim-k; l++)
	fprintf(f, "  %e %e",
		creal(MEL(dim, mat, l, l+k)), cimag(MEL(dim, mat, l, l+k)));
    fprintf(f, "\n");
  }
  fclose(f);

  /*
   * greens_function.out
   */
  f = fopen(gffn, "w");
  if (!f) {
    printf("Error opening %s.\n", gffn);
    exit(999);
  }
  fprintf(f, "# dist  qx qy");
  dim = kernel->get_dimension();
  i = 4;
  for (int k = 0; k < dim; k++)
    for (int l = 0; l < dim-k; l++) {
      fprintf(f, "  %i:mat%i%i", i, l+1, l+k+1);
      i += 2;
    }
  fprintf(f, "\n");
  for (i = 0; i < npts; i++) {
    kernel->get_stiffness_matrix(2*M_PI*ptx[i], 2*M_PI*pty[i], mat);
    GaussJordan(dim, mat, &error);
    fprintf(f, "%e  %e %e", ptd[i], ptx[i], pty[i]);
    for (int k = 0; k < dim; k++)
      for (int l = 0; l < dim-k; l++)
	fprintf(f, "  %e %e",
		creal(MEL(dim, mat, l, l+k)), cimag(MEL(dim, mat, l, l+k)));
    fprintf(f, "\n");
  }
  fclose(f);


  f = fopen("specpoints.out", "w");
  fprintf(f, "# symbol dist qx qy\n");
  for (i = 0; i < nspec; i++) {
    //    fprintf(f, "%c  %e  %e %e\n", specsym[i], specd[i], specx[i], specy[i]);
    fprintf(f, "%s  %e  %e %e\n", &specsym[i*SYM_LEN], specd[i], specx[i],
	    specy[i]);
  }
  fclose(f);

  free(specsym);
  free(specx);
  free(specy);
  free(specd);
  free(ptx);
  free(pty);
  free(ptd);
}
