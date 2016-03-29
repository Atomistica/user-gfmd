#ifndef __GEOMETRY_H
#define __GEOMETRY_H

double rand1();

double clockwise_angle(double x1, double y1, double x2, double y2);

bool lines_cross(double A1x, double A1y, double A2x, double A2y,
                 double B1x, double B1y, double B2x, double B2y);

bool tr_is_inside(double Bx, double By, double A1x, double A1y,
                  double A2x, double A2y, double A3x, double A3y);
                  
double tr_area(double A1x, double A1y, double A2x, double A2y, double A3x,
               double A3y);

bool tr_overlap(double Bx, double By, double Cx, double Cy,
                double A1x, double A1y, double A2x, double A2y);

#endif
