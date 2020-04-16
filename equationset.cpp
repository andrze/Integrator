#include "equationset.h"
#include "physics_cubic.h"
#include <cmath>
#include <iostream>

EquationSet::EquationSet(bool cubic, double dimension)
{
    this->cubic = cubic;
    this->d = dimension;
    d_factor = std::pow(2, -1-dimension)*std::pow(M_PI,-d/2)/std::tgamma(dimension/2); //d_factor = v_d from Delamotte's convention
}

RealVector EquationSet::evaluate(RealVector point){
    //if(cubic){

    double eta = find_eta(point, d, d_factor);
    std::vector<double> coords;
    coords.push_back(K_der_cubic(point, eta, d, d_factor));
    coords.push_back(Z_der_cubic(point, eta));
    coords.push_back(Y_der_cubic(point, eta, d, d_factor));
    coords.push_back(u20_der_cubic(point, eta, d, d_factor));
    coords.push_back(u01_der_cubic(point, eta, d, d_factor));
    coords.push_back(u30_der_cubic(point, eta, d, d_factor));
    coords.push_back(u11_der_cubic(point, eta, d, d_factor));
    coords.push_back(u40_der_cubic(point, eta, d, d_factor));
    coords.push_back(u21_der_cubic(point, eta, d, d_factor));
    coords.push_back(u02_der_cubic(point, eta, d, d_factor));

    return RealVector(coords);

    //}
}
