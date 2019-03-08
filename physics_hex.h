#ifndef PHYSICS_HEX_H
#define PHYSICS_HEX_H
#include "realvector.h"

double A_der_hex(RealVector x, double L, double zp_der, double d, double d_factor);
double Ms_der_hex(RealVector x, double L, double zp_der, double d, double d_factor);
double Mp_der_hex(RealVector x, double L, double zp_der, double d, double d_factor);
double Zs_der_hex(RealVector x, double L, double zp_der, double d, double d_factor);
double Zp_der_hex(RealVector x, double L, double d, double d_factor);


#endif // PHYSICS_HEX_H
