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
    EquationSet(std::vector<std::function<double (RealVector, double, double)> > equations,
                std::vector<std::function<double (double, double)> > scale);

    std::vector<std::function<double (RealVector, double, double ) > > equations;
    std::vector<std::function<double (double, double)> >  scale;

    RealVector evaluate(RealVector point, double d);
};

#endif // REALVECTOR_H
