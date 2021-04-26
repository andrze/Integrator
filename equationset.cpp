#include "equationset.h"
#include "physics.h"
#include <cmath>
#include <iostream>

EquationSet::EquationSet(double dimension,double n)
{
    this->d = dimension;
    this->n = n;
    d_factor = std::pow(2, -1-dimension)*std::pow(M_PI,-d/2)/std::tgamma(dimension/2); //d_factor = v_d from Delamotte's convention
}

RealVector EquationSet::evaluate(RealVector point){
    //if(cubic){

    double eta = find_eta(point, d, n);
    std::vector<double> coords;
    coords.push_back(K_der(point, eta, d, n));
    coords.push_back(u20_der(point, eta, d, n));
    coords.push_back(Z_der(point, eta));

    return RealVector(coords);

    //}
}
