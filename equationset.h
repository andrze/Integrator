#ifndef EQUATIONSET_H
#define EQUATIONSET_H
#include <vector>
#include "realvector.h"
#include "functional"


class EquationSet
{
public:
    EquationSet(bool cubic=true);

    bool cubic;

    RealVector evaluate(RealVector point, double L);
};
#endif // EQUATIONSET_H
