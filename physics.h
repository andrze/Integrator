#ifndef PHYSICS_CUBE_H
#define PHYSICS_CUBE_H
#include "realvector.h"

double K_der(RealVector x, double eta, double d, double n);

double Z_der(RealVector point, double eta);
double find_eta(RealVector x, double d, double d_factor);

double u20_der(RealVector x, double eta, double d, double n);



#endif // PHYSICS_CUBE_H
