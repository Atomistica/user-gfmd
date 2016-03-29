#ifndef __LI_BERGER_H
#define __LI_BERGER_H

namespace LiBerger {

void tr_uniform_PQ_c(double G, double nu,
                     double P, double Qx, double Qy, double Bx, double By,
		     double A2x, double A2y, double A3x, double A3y,
		     double &ux, double &uy, double &uz);

void tr_uniform_PQ(double G, double nu,
		   double P, double Qx, double Qy, double A1x, double A1y,
		   double A2x, double A2y, double A3x, double A3y,
		   double Bx, double By, double &ux, double &uy, double &uz);
		           
void tr_linear_PQ_c(double G, double nu,
                    double P, double Qx, double Qy, double Bx, double By,
                    double A2x, double A2y, double A3x, double A3y,
                    double &ux, double &uy, double &uz);

void tr_linear_PQ(double G, double nu,
                  double P, double Qx, double Qy, double A1x, double A1y,
                  double A2x, double A2y, double A3x, double A3y,
                  double Bx, double By, double &ux, double &uy, double &uz);

void sq_uniform_PQ(double G, double nu,
		   double P, double Qx, double Qy, double a, double b, 
		   double Bx, double By, double &ux, double &uy, double &uz);

}

#endif
