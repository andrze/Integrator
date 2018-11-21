#ifndef PHYSICS_HEX_H
#define PHYSICS_HEX_H
#include "realvector.h"

double A_der_hex(RealVector x, double L, double zp_der);
double Ms_der_hex(RealVector x, double L, double zp_der);
double Mp_der_hex(RealVector x, double L, double zp_der);
double Zs_der_hex(RealVector x, double L, double zp_der);
double Zp_der_hex(RealVector x, double L);


#endif // PHYSICS_HEX_H
