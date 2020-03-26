#ifndef EQUATIONSET_H
#define EQUATIONSET_H
#include <vector>
#include "realvector.h"


class EquationSet
{
public:
    EquationSet(bool cubic=true, double dimension=2);

    bool cubic;
    double d;
    double d_factor;

    RealVector evaluate(RealVector point);
};
#endif // EQUATIONSET_H
