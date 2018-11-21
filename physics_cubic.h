#ifndef PHYSICS_CUBE_H
#define PHYSICS_CUBE_H
#include "realvector.h"

double A_der_cubic(RealVector x, double L, double zp_der);
double Ms_der_cubic(RealVector x, double L, double zp_der);
double Mp_der_cubic(RealVector x, double L, double zp_der);
double Zs_der_cubic(RealVector x, double L, double zp_der);
double Zp_der_cubic(RealVector x, double L);


#endif // PHYSICS_CUBE_H
