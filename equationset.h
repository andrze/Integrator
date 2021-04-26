#ifndef EQUATIONSET_H
#define EQUATIONSET_H
#include <vector>
#include "realvector.h"


class EquationSet
{
public:
    EquationSet(double dimension=2, double n=2);

    double d;
    double d_factor;
    double n;

    RealVector evaluate(RealVector point);
};
#endif // EQUATIONSET_H
