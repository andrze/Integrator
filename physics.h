#ifndef PHYSICS_H
#define PHYSICS_H
#include "realvector.h"

double A_der(RealVector x, double L, double zp_der);
double Ms_der(RealVector x, double L, double zp_der);
double Mp_der(RealVector x, double L, double zp_der);
double Zs_der(RealVector x, double L, double zp_der);
double Zp_der(RealVector x, double L);

#endif // PHYSICS_H
