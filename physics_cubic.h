#ifndef PHYSICS_CUBE_H
#define PHYSICS_CUBE_H
#include "realvector.h"

double K_der_cubic(RealVector x, double eta, double d, double d_factor);

double Z_der_cubic(RealVector point, double eta);
double find_eta(RealVector x, double d, double d_factor);
double Y_der_cubic(RealVector x, double eta, double d, double d_factor);

double u20_der_cubic(RealVector x, double eta, double d, double d_factor);
double u01_der_cubic(RealVector x, double eta, double d, double d_factor);

double u30_der_cubic(RealVector x, double eta, double d, double d_factor);
double u11_der_cubic(RealVector x, double eta, double d, double d_factor);

double u40_der_cubic(RealVector x, double eta, double d, double d_factor);
double u21_der_cubic(RealVector x, double eta, double d, double d_factor);
double u02_der_cubic(RealVector x, double eta, double d, double d_factor);





#endif // PHYSICS_CUBE_H
