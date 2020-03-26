#include "equationset.h"
#include "physics_cubic.h"
#include <cmath>
#include <iostream>

EquationSet::EquationSet(bool cubic, double dimension)
{
    this->cubic = cubic;
    this->d = dimension;
    d_factor = std::pow(2, -1-dimension)*std::pow(M_PI,-d/2)/std::tgamma(dimension/2);
}

RealVector EquationSet::evaluate(RealVector point){
    //if(cubic){

    double eta = find_eta(point, d, d_factor);
    std::vector<double> coords;
    coords.push_back(K_der_cubic(point, eta, d, d_factor));
    coords.push_back(u_der_cubic(point, eta, d, d_factor));
    coords.push_back(l_der_cubic(point, eta, d, d_factor));
    coords.push_back(v_der_cubic(point, eta, d, d_factor));
    coords.push_back(Y_der_cubic(point, eta, d, d_factor));
    coords.push_back(Z_der_cubic(point, eta));
    coords.push_back(J_der_cubic(point, eta, d, d_factor));

    return RealVector(coords);

    //}
}
