#ifndef PHYSICS_CUBE_H
#define PHYSICS_CUBE_H
#include "realvector.h"

double A_der_cubic(RealVector x, double L, double zp_der, double d, double d_factor);
double Ms_der_cubic(RealVector x, double L, double zp_der, double d, double d_factor);
double Mp_der_cubic(RealVector x, double L, double zp_der, double d, double d_factor);
double Zs_der_cubic(RealVector x, double L, double zp_der, double d, double d_factor);
double Zp_der_cubic(RealVector x, double L, double d, double d_factor);


#endif // PHYSICS_CUBE_H
