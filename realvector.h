#ifndef REALVECTOR_H
#define REALVECTOR_H

#include <vector>
#include <ostream>
#include <functional>

class RealVector
{
public:
    RealVector();
    RealVector(std::vector<double> coords);


    std::vector<double> coords;

    RealVector& operator += (RealVector rhs);
    RealVector& operator *= (double rhs);
    RealVector& operator -= (RealVector rhs);
    RealVector& operator /= (double rhs);
    double& operator [] (size_t i);
};

RealVector operator + (RealVector lhs, RealVector rhs);
RealVector operator - (RealVector lhs, RealVector rhs);
RealVector operator * (RealVector lhs, double rhs);
RealVector operator * (double lhs, RealVector rhs);
RealVector operator / (RealVector lhs, double rhs);
double operator * (RealVector lhs, RealVector rhs);
std::ostream& operator << (std::ostream& out, RealVector v);

RealVector filter(RealVector v, std::vector<bool>);

struct EquationSet
{
    EquationSet();
    EquationSet(std::vector<std::function<double(RealVector ) > > equations,
                std::vector<bool> differential);

    std::vector<std::function<double(RealVector ) > > equations;
    std::vector<bool> differential;

    RealVector evaluate(RealVector point);
};

#endif // REALVECTOR_H
