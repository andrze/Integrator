#ifndef PHYSICS_CUBE_H
#define PHYSICS_CUBE_H
#include "realvector.h"

double K_der_cubic(RealVector x, double eta, double d, double d_factor);
double u_der_cubic(RealVector x, double eta, double d, double d_factor);
double l_der_cubic(RealVector x, double eta, double d, double d_factor);
double v_der_cubic(RealVector x, double eta, double d, double d_factor);
double Y_der_cubic(RealVector x, double eta, double d, double d_factor);
double Z_der_cubic(RealVector point, double eta);
double J_der_cubic(RealVector x, double eta, double d, double d_factor);

double find_eta(RealVector x, double d, double d_factor);


#endif // PHYSICS_CUBE_H
