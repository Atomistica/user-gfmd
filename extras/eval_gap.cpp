/*
 * Standalone code, count number of atoms that penetrate a contact map.
 */
#include <cmath>

#include <stdio.h>

#include <netcdfcpp.h>

#include "error.h"
#include "memory.h"
#include "table2d.h"

#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))

using namespace std;
using namespace LAMMPS_NS;

enum { AVERAGE, BOTTOM, TOP };

void read_map(char *fn, int z0_type, double z0, long &nx, long &ny, double **&map_data,
	      Memory *memory, Error *error)
{
  FILE *f = fopen(fn, "r");

  if (!f) {
    char errstr[1024];
    sprintf(errstr, "Could not open contact map file '%s'.", fn);
    error->all(errstr);
  }

  // count number of columns
  nx = 0;
  char c = fgetc(f), lastc = ' ';
  while (!feof(f) && c != '\n') {
    if ((isalnum(c) || c == '-' || c == 'e') && isblank(lastc))
      nx++;
    lastc = c;
    c = fgetc(f);
  }
  if (feof(f))
    error->all("fix contact/map: Could not read map file.");

  if (nx <= 0)
    error->all("fix contact/map: Something wrong: nx <= 0");


  /*
   * for now assume a square map
   */
  ny = nx;

  rewind(f);

  memory->create(map_data, nx, ny, "map");
  double za = 0.0, maxz = map_data[0][0], minz = map_data[0][0];
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      double tmp;
      if (fscanf(f, "%lf", &tmp) != 1)
	error->all("fix contact/map: Error read map file.");

      map_data[i][j] = tmp;

      za += tmp;

      maxz = MAX(maxz, tmp);
      minz = MIN(minz, tmp);
    }
  }

  fclose(f);

  // remove average height and add zero height
  switch (z0_type) {
  case BOTTOM:
    za = z0 - minz;
    break;
  case TOP:
    za = z0 - maxz;
    break;
  default:
    za = z0 - za/(nx*ny);
    break;
  }

  for (int i = 0; i < nx*ny; i++) {
    map_data[0][i] += za;
  }
}


void write_map(char *fn, long nx, long ny, double **map_data, Error *error)
{
  FILE *f = fopen(fn, "w");

  if (!f) {
    char errstr[1024];
    sprintf(errstr, "Could not open contact map file '%s'.", fn);
    error->all(errstr);
  }

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      fprintf(f, "%lf ", map_data[i][j]);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}


double min_dist(int w, int x, int y, long nx, long ny, double **gap_data)
{
  double d = w;
  int w2 = w/2;

  //  printf("%i %i %i %i %i\n", w, x, y, nx, ny);

  for (int i = -w2; i < w2; i++) {
    for (int j = -w2; j < w2; j++) {
      long ii = x+i;
      long jj = y+j;

      while (ii <  0)   ii += nx;
      while (ii >= nx)  ii -= nx;
      while (jj <  0)   jj += ny;
      while (jj >= ny)  jj -= ny;

      //      printf("%i %i\n", ii, jj);
      if (gap_data[ii][jj] < 1.0) {
	d = MIN(d, sqrt(i*i + j*j));
      }
    }
  }

  return d;
}


int main(int narg, char **argv)
{
  Error *error = new Error();
  Memory *memory = new Memory(error);

  double z0, ming, maxg, maxr, dg;
  int ng, nr;
  char *endptr;
  bool write_gap_map = false;

  if (narg < 7) {
    printf("Syntax: eval_gap <trajectory> <map> <z0> <ming> <maxg> <ng> <maxr> <nr> [gap-map]"
	   "[<h1>] [<h2>] [<h3>] ...\n");
    return 999;
  }

  int z0_type;
  int iarg = 3;
  z0_type = AVERAGE;
  if (!strcmp(argv[iarg], "bot")) {
    z0_type = BOTTOM;
    iarg++;
  }
  else if (!strcmp(argv[iarg], "top")) {
    z0_type = TOP;
    iarg++;
  }

  z0 = strtod(argv[iarg], &endptr);
  if (endptr == argv[iarg])
    error->all("Could not convert z0 argument to double.");
  iarg++;
  ming = strtod(argv[iarg], &endptr);
  if (endptr == argv[iarg])
    error->all("Could not convert ming argument to double.");
  iarg++;
  maxg = strtod(argv[iarg], &endptr);
  if (endptr == argv[iarg])
    error->all("Could not convert maxg argument to double.");
  iarg++;
  ng = strtol(argv[iarg], &endptr, 10);
  if (endptr == argv[iarg])
    error->all("Could not convert ng argument to double.");
  iarg++;
  maxr = strtod(argv[iarg], &endptr);
  if (endptr == argv[iarg])
    error->all("Could not convert maxr argument to double.");
  iarg++;
  nr = strtol(argv[iarg], &endptr, 10);
  if (endptr == argv[iarg])
    error->all("Could not convert nr argument to double.");
  iarg++;

  if (iarg < narg) {
    if (!strcmp(argv[iarg], "gap-map")) {
      write_gap_map = true;
      iarg++;
    }
  }

  dg = (maxg-ming)/ng;

  double *h = NULL;
  int *hn = NULL;
  int nh = narg-iarg;
  if (nh > 0) {
    h = new double[nh];
    hn = new int[nh];
    for (int i = iarg; i < narg; i++) {
      h[i-iarg] = strtod(argv[i], &endptr);
      if (endptr == argv[i])
	error->all("Could not convert h argument to double.");
      printf("%i %i %f\n", i-iarg, nh, h[i-iarg]);
    }
  }

  long nx, ny;
  double **map_data;
  read_map(argv[2], z0_type, z0, nx, ny, map_data, memory, error);
  Table2D *map = new Table2D(nx, ny, map_data, true, error, memory);
  memory->destroy(map_data);

  NcFile *nc = new NcFile(argv[1]);

  int *hist = new int[ng];
  double *avgd = new double[nr];
  double *avgg = new double[nr];
  double *devg = new double[nr];
  int *numg = new int[nr];

  int nframes = nc->get_dim("frame")->size();
  int natoms = nc->get_dim("atom")->size();

  int nelx = int(sqrt(natoms));
  int nely = nelx;

  if (nelx*nely != natoms) {
    error->all("Atoms not on a square grid.");
  }

  int *ix = NULL, *iy = NULL;
  memory->create(map_data, nelx, nely, "map");

  NcVar *x_var = nc->get_var("coordinates");
  double **x;
  memory->create(x, natoms, 3, "x");
  double *x_ptr = x[0];

  NcVar *l_var = nc->get_var("cell_lengths");

  FILE *hf = NULL;
  if (h) {
    hf = fopen("area_from_gap.out", "w");
    fprintf(hf, "#");
    for (int j = 0; j < nh; j++) {
      fprintf(hf, " %f", h[j]);
    }
    fprintf(hf, "\n");
  }

  if (write_gap_map) {
    printf("Writing gap map for each frame.\n");
  }

  for (int frame = 0; frame < nframes; frame++) {
    memset(hist, 0, ng*sizeof(int));
    memset(avgd, 0, nr*sizeof(double));
    memset(avgg, 0, nr*sizeof(double));
    memset(devg, 0, nr*sizeof(double));
    memset(numg, 0, nr*sizeof(int));

    if (!x_var->set_cur(frame))
      error->all("set_cur failed.");
    if (!x_var->get(x_ptr, 1, natoms, 3))
      error->all("get failed.");

    double l[3];
    if (!l_var->set_cur(frame))
      error->all("set_cur failed.");
    if (!l_var->get(l, 1, 3))
      error->all("get failed.");

    if (!ix) {
      memory->create(ix, natoms, "ix");
      memory->create(iy, natoms, "iy");

      memset(*map_data, 0, nelx*nely*sizeof(double));

      for (int i = 0; i < natoms; i++) {
	ix[i] = int(x[i][0]*nelx/l[0]) % nelx;
	iy[i] = int(x[i][1]*nely/l[1]) % nely;

	if (map_data[ix[i]][iy[i]] > 0.0) {
	  char errstr[1024];
	  sprintf(errstr, "Atom at map entry %i x %i exists twice.", ix[i], iy[i]);
	  error->all(errstr);
	}

	map_data[ix[i]][iy[i]] = 1.0;
      }

      for (int i = 0; i < nelx; i++) {
	for (int j = 0; j < nely; j++) {
	  if (map_data[i][j] < 1.0) {
	    char errstr[1024];
	    sprintf(errstr, "Could not find atom for map entry %i x %i.", i, j);
	    error->all(errstr);
	  }
	}
      }
    }

    double rx = x[0][0]*nx/l[0];
    double ry = x[0][1]*ny/l[1];
    double z, dz_dx, dz_dy;
    map->eval(rx, ry, z, dz_dx, dz_dy);
    z -= x[0][2];

    double cur_ming = z, cur_maxg = z;

    if (h) {
      memset(hn, 0, nh*sizeof(int));
    }

    int n = 0;
    for (int i = 0; i < natoms; i++) {
      rx = x[i][0]*nx/l[0];
      ry = x[i][1]*ny/l[1];

      map->eval(rx, ry, z, dz_dx, dz_dy);

      z -= x[i][2];

      map_data[ix[i]][iy[i]] = z;

      if (abs(z) < 1.0) {
	n++;
      }

      if (h) {
	for (int j = 0; j < nh; j++) {
	  if (abs(z) < h[j])
	    hn[j]++;
	}
      }

      int j = int((z-ming)/dg);
      if (j >= 0 && j < ng)
	hist[j]++;
      cur_ming = MIN(cur_ming, z);
      cur_maxg = MAX(cur_maxg, z);
    }

    for (int i = 0; i < natoms; i++) {
      double g = map_data[ix[i]][iy[i]];
      if (g >= 1.0) {
	double d = min_dist(int(maxg)+1, ix[i], iy[i], nelx, nely, map_data);

	int j = int(d*nr/maxr);

	if (j < nr) {
	  avgd[j] += d;
	  avgg[j] += g;
	  devg[j] += g*g;
	  numg[j]++;
	}
      }
    }

    for (int j = 0; j < nr; j++) {
      if (numg[j] > 0) {
	avgd[j] /= numg[j];
	avgg[j] /= numg[j];
	devg[j]  = sqrt(devg[j]/numg[j] - avgg[j]*avgg[j]);
      }
    }

    printf("%i  (%f %f) - %f\n", frame, cur_ming, cur_maxg, n/(l[0]*l[1]));

    if (h) {
      fprintf(hf, "%i", frame);
      for (int j = 0; j < nh; j++) {
	fprintf(hf, " %f", ((double) hn[j])/(nelx*nely));
      }
      fprintf(hf, "\n");
    }

    char fn[1024];
    sprintf(fn, "%i.gap.out", frame);
    FILE *f = fopen(fn, "w");
    double A = n/(l[0]*l[1]);
    fprintf(f, "# A = %f, A0 = %f, natoms = %i\n", A, l[0]*l[1], natoms);
    for (int i = 0; i < ng; i++) {
      fprintf(f, "%e %e %e %i\n", ming+(i+0.5)*dg,
	      ((double)hist[i])/natoms/dg, ((double)hist[i])/(natoms*A)/dg, hist[i]);
    }
    fclose(f);

    sprintf(fn, "%i.perimeter.out", frame);
    f = fopen(fn, "w");
    fprintf(f, "# A = %f\n", A);
    fprintf(f, "# 1:distance 2:gap 3:stddev-gap 4:numg 5:nP 6:nP/A0 7:nP/Arep\n");
    int nP = 0;
    for (int i = 0; i < nr; i++) {
      if (numg[i] > 0) {
      //      fprintf(f, "%e %e %e %i\n", (i+0.5)*maxr/nr, avgg[i], devg[i], numg[i]);
	nP += numg[i];
	fprintf(f, "%e %e %e %i %i %e %e\n", avgd[i], avgg[i], devg[i], numg[i], nP, nP/(l[0]*l[1]), float(nP)/float(n));
      }
    }
    fclose(f);

    if (write_gap_map) {
      sprintf(fn, "%i.gap_map.out", frame);
      write_map(fn, nelx, nely, map_data, error);
    }
  }

  if (h) {
    fclose(hf);
    delete [] h;
    delete [] hn;
  }

  if (ix) {
    memory->destroy(ix);
    memory->destroy(iy);
  }
  memory->destroy(map_data);
  memory->destroy(x);

  delete hist;
  delete avgg;
  delete devg;
  delete numg;

  delete nc;
  //  delete map;
  delete memory;
  delete error;
}
