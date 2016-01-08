/*
 * Boussinesq-Cerruti solution for uniform pressure on a triangular area.
 * Formulae and notation is from: Li, Berger, J. Elast. 63, 137 (2001)
 *
 * Sign convention: negative P is pushing down, positive P is pushing up.
 * Similarly: negative uz is displacement into the surface.
 * I.e. signs are reverse to what is usually assumed (Johnson; Li, Berger)
 */

#include <assert.h> 
#include <math.h>
#include <stdio.h>

#include "geometry.h"

#include "li_berger.h"

namespace LiBerger {

/*
 * These are the fundamental integrals
 */
 
double _Ii100(double phi)
{
    return log(fabs(tan(0.5*phi)));
}

double Ii100(double n1, double eps1, double eps2)
{
    return n1*(_Ii100(eps1) - _Ii100(eps2));
}


double _Ii101(double delta1, double phi)
{
    return sin(delta1)/sin(phi) + cos(delta1)*log(fabs(tan(0.5*phi)));
}

double Ii101(double n1, double delta1, double eps1, double eps2)
{
    return 0.5*n1*n1*(_Ii101(delta1, eps1) - _Ii101(delta1, eps2));
}


double Ii201(double n1, double delta1, double eps1, double eps2)
{
    double sdi = sin(delta1);
    double cdi = cos(delta1);

    double t1 = sdi*log(fabs(sin(eps1))) - eps1*cdi;
    double t2 = sdi*log(fabs(sin(eps2))) - eps2*cdi;

    return -n1*(t1-t2);
}


double Ii210(double n1, double delta1, double eps1, double eps2)
{
    double sdi = sin(delta1);
    double cdi = cos(delta1);

    double t1 = cdi*log(fabs(sin(eps1))) + eps1*sdi;
    double t2 = cdi*log(fabs(sin(eps2))) + eps2*sdi;

    return -n1*(t1-t2);
}


double _Ii311(double delta1, double phi)
{
    return sin(2*delta1)*log(fabs(tan(0.5*phi))) +
           2*sin(2*delta1)*cos(phi) -
           2*cos(2*delta1)*sin(phi);
}

double Ii311(double n1, double delta1, double eps1, double eps2)
{
    return 0.5*n1*(_Ii311(delta1, eps1) - _Ii311(delta1, eps2));
}


/* Note: Ii302 is not in Li, Berger */
double _Ii302(double delta1, double phi)
{
    double sdi = sin(delta1);

    return sdi*sdi*log(fabs(tan(0.5*phi))) - cos(2*delta1-phi);
}

double Ii302(double n1, double delta1, double eps1, double eps2)
{
    return n1*(_Ii302(delta1, eps1) - _Ii302(delta1, eps2));
}


/* Note: Ii320 is identical expression but different formulation than
   Li, Berger */
double _Ii320(double delta1, double phi)
{
    double cdi = cos(delta1);

    return cdi*cdi*log(fabs(tan(0.5*phi))) + cos(2*delta1-phi);
}

double Ii320(double n1, double delta1, double eps1, double eps2)
{
    return n1*(_Ii320(delta1, eps1) - _Ii320(delta1, eps2));
}


double _Ii331(double delta1, double phi)
{
    return -0.5*(sin(2*delta1)+0.5*sin(4*delta1))
           *(cos(phi)/(sin(phi)*sin(phi))+asinh(1.0/tan(phi)))
           -(2*sin(2*delta1)+4*sin(4*delta1))*log(fabs(0.5*phi))
           +(2*cos(2*delta1)+2*cos(4*delta1))/sin(phi)
           -4*sin(4*delta1)*cos(phi)
           +4*cos(4*delta1)*sin(phi);
}

double Ii331(double n1, double delta1, double eps1, double eps2)
{
    return n1*n1*n1/12*(_Ii311(delta1, eps1) - _Ii311(delta1, eps2));
}


/*!
 * Compute displacements at point Bx, By due to uniform pressure P
 * and tractions Qx, Qy on a triangle bounded by Bx,By, A2x,A2y and A3x,A3y.
 * The point where the displacement is computed sits at the corner of the
 * triangle.
 */
void tr_uniform_PQ_c(double G, double nu,
                     double P, double Qx, double Qy, double Bx, double By,
		     double A2x, double A2y, double A3x, double A3y,
		     double &ux, double &uy, double &uz)
{
    /* Area from cross product, area = 1/2 | (A2-B) x (A3-B) | */
    double area = 0.5*fabs( (A2x-Bx)*(A3y-By) - (A3x-Bx)*(A2y-By) );
    
    if (area <= 0.0)
        return;
    
    /* Length of the sides of the triangle */
    double lBA2 = sqrt( (Bx-A2x)*(Bx-A2x) + (By-A2y)*(By-A2y) );
    double lBA3 = sqrt( (Bx-A3x)*(Bx-A3x) + (By-A3y)*(By-A3y) );
    double lA2A3 = sqrt( (A2x-A3x)*(A2x-A3x) + (A2y-A3y)*(A2y-A3y) );
    
    /* height of the triangle */
    double n1 = 2*area/lA2A3;
    /* angle with respect to x-axis */
    double delta1 = -clockwise_angle(1.0,0.0, A3x-A2x,A3y-A2y);
    /* epsilon angles, see Figs. 2 and 3 of Li & Berger */
    double eps1 = clockwise_angle(A2x-A3x,A2y-A3y, Bx-A2x,By-A2y);
    double eps2 = clockwise_angle(A2x-A3x,A2y-A3y, Bx-A3x,By-A3y);
    
    /*
     * Evaluate integrals
     */
    double v_Ii100 = Ii100(n1, eps1, eps2);
    double v_Ii201 = Ii201(n1, delta1, eps1, eps2);
    double v_Ii210 = Ii210(n1, delta1, eps1, eps2);
    double v_Ii311 = Ii311(n1, delta1, eps1, eps2);
    double v_Ii302 = Ii302(n1, delta1, eps1, eps2);
    double v_Ii320 = Ii320(n1, delta1, eps1, eps2);
    
    /*
     * Displacements due to normal load P
     */
    ux += P*(1-2*nu)/(4*M_PI*G) * v_Ii210;
    uy += P*(1-2*nu)/(4*M_PI*G) * v_Ii201;
    uz += P*(1-nu)/(2*M_PI*G) * v_Ii100;

    /*
     * Displacements due to traction Qx
     */
    ux +=  Qx*( (1-nu)/(2*M_PI*G)*v_Ii100 + nu/(2*M_PI*G)*v_Ii320 );
    uy +=  Qx*nu/(2*M_PI*G)*v_Ii311;
    uz += -Qx*(1-2*nu)/(4*M_PI*G)*v_Ii210;
  
    /*
     * Displacements due to traction Qy
     */
    ux +=  Qy*nu/(2*M_PI*G)*v_Ii311;
    uy +=  Qy*( (1-nu)/(2*M_PI*G)*v_Ii100 + nu/(2*M_PI*G)*v_Ii302 );
    uz += -Qy*(1-2*nu)/(4*M_PI*G)*v_Ii201;
}


/*!
 * Compute displacements at point Bx, By due to uniform pressure P
 * and tractions Qx, Qy on a triangle bounded by A1x,A1y, A2x,A2y and A3x,A3y.
 * See Fig. 2 of Li & Berger. Appropriate signs need to be attached if point is
 * outside.
 */
void tr_uniform_PQ(double G, double nu,
		   double P, double Qx, double Qy, double A1x, double A1y,
		   double A2x, double A2y, double A3x, double A3y,
		   double Bx, double By, double &ux, double &uy, double &uz)
{
    double area = tr_area(A1x,A1y, A2x,A2y, A3x,A3y);
    double area_tot = 0.0;
    
    area = tr_area(Bx,By, A1x,A1y, A2x,A2y);
    if (tr_overlap(Bx,By, A3x,A3y, A1x,A1y, A2x,A2y)) {
        tr_uniform_PQ_c(G,nu, P,Qx,Qy, Bx,By, A1x,A1y, A2x,A2y, ux,uy,uz);
        area_tot += area;
    }
    else {
        tr_uniform_PQ_c(G,nu, -P,-Qx,-Qy, Bx,By, A1x,A1y, A2x,A2y, ux,uy,uz);
        area_tot -= area;
    }

    area = tr_area(Bx,By, A2x,A2y, A3x,A3y);
    if (tr_overlap(Bx,By, A1x,A1y, A2x,A2y, A3x,A3y)) {
        tr_uniform_PQ_c(G,nu, P,Qx,Qy, Bx,By, A2x,A2y, A3x,A3y, ux,uy,uz);
        area_tot += area;
    }
    else {
        tr_uniform_PQ_c(G,nu, -P,-Qx,-Qy, Bx,By, A2x,A2y, A3x,A3y, ux,uy,uz);
        area_tot -= area;
    }

    area = tr_area(Bx,By, A3x,A3y, A1x,A1y);    
    if (tr_overlap(Bx,By, A2x,A2y, A3x,A3y, A1x,A1y)) {
        tr_uniform_PQ_c(G,nu, P,Qx,Qy, Bx,By, A3x,A3y, A1x,A1y, ux,uy,uz);
        area_tot += area;
    }
    else {
        tr_uniform_PQ_c(G,nu, -P,-Qx,-Qy, Bx,By, A3x,A3y, A1x,A1y, ux,uy,uz);
        area_tot -= area;
    }

    assert(fabs(tr_area(A1x,A1y, A2x,A2y, A3x,A3y)-area_tot) < 1e-9);
}


/*!
 * Compute displacements at point Bx, By due to uniform pressure P
 * and tractions Qx, Qy on a rectangular area centered at the origin with
 * sides of length 2a,2b.
 */
void sq_uniform_PQ(double G, double nu,
		   double P, double Qx, double Qy, double a, double b, 
		   double Bx, double By, double &ux, double &uy, double &uz)
{
    tr_uniform_PQ(G,nu, P,Qx,Qy, -a,-b, a,-b, a,b, Bx,By, ux,uy,uz);
    tr_uniform_PQ(G,nu, P,Qx,Qy, -a,-b, -a,b, a,b, Bx,By, ux,uy,uz);
}

}
