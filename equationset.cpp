#include "equationset.h"
#include "physics_cubic.h"
#include "physics_hex.h"
#include <cmath>

EquationSet::EquationSet(bool cubic, double dimension)
{
    this->cubic = cubic;
    this->d = dimension;
    d_factor = 2/d*std::pow(M_PI, dimension/2)*std::pow(2*M_PI,-d)/std::tgamma(dimension/2);
}

RealVector EquationSet::evaluate(RealVector point, double L){
    if(cubic){
        double zp_der = Zp_der_cubic(point, L, d, d_factor);
        std::vector<double> coords;
        coords.push_back(A_der_cubic(point, L, zp_der, d, d_factor));
        coords.push_back(Ms_der_cubic(point, L, zp_der, d, d_factor));
        coords.push_back(Mp_der_cubic(point, L, zp_der, d, d_factor));
        coords.push_back(Zs_der_cubic(point, L, zp_der, d, d_factor));
        coords.push_back(zp_der);
        coords.push_back(0.);

        return RealVector(coords);

    } else {
    	double zp_der = Zp_der_hex(point, L, d, d_factor);
		std::vector<double> coords;
		coords.push_back(A_der_hex(point, L, zp_der, d, d_factor));
		coords.push_back(Ms_der_hex(point, L, zp_der, d, d_factor));
		coords.push_back(Mp_der_hex(point, L, zp_der, d, d_factor));
		coords.push_back(Zs_der_hex(point, L, zp_der, d, d_factor));
		coords.push_back(zp_der);
		coords.push_back(0.);

		return RealVector(coords);
    }
}
