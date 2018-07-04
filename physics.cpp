#include <cmath>
#include <vector>
#include "realvector.h"
#include "physics.h"

//extern const double d = 3., gamma3 = gamma(3);

double I(double d, int n, double w, double eta=0){

    return (1-eta/(d+2))/(std::pow(1+w,n));
}

double mL(RealVector x){
    return 2*x[0]*x[1];
}

double mT(RealVector x){
    return 18*x[0]*x[0]*x[2];
}

double k_der(RealVector x, double eta, double d){
    double dimensional = (2-d-eta)*x[0];
    double nondimensional;

    if( (x[0]*x[2]==0) && (x[1] == 0)){
        nondimensional = 3*I(d, 2, mL(x),eta) + I(d,2,mT(x),eta);
    } else {
        nondimensional = 3*I(d, 2, mL(x),eta) + (1+36*x[0]*x[2]/x[1])*I(d,2,mT(x),eta);
    }

    return dimensional + nondimensional;
}

double k_scale(double eta, double d){
    return d-2+eta;
}

double u_der(RealVector x, double eta, double d){
    double dimensional = (d-4-2*eta)*x[1];
    double isotropic = 18*x[1]*x[1]*I(d,3,mL(x),eta);
    double anisotropic = 2*(x[1]+36*x[0]*x[2])*(x[1]+36*x[0]*x[2]) * I(d,3,mT(x),eta) -
            36*x[2] * I(d,2,mT(x),eta);
    return dimensional + isotropic + anisotropic;
}

double u_scale(double eta, double d){
    return d-4+2*eta;
}

double l_der(RealVector x, double eta, double d){
    if(x[2] == 0){
        return 0;
    }
    double dimensional = (2*d-6-3*eta)*x[2];
    double nondimensional = 15*x[2]*(x[1]+6*x[0]*x[2])/(x[0]*(x[1]-9*x[0]*x[2])) *
            (I(d,2,mT(x),eta) -  I(d,2,mL(x),eta));
    return dimensional + nondimensional;
}

double l_scale(double eta, double d){
    return 2*d-6+3*eta;
}

double eta(RealVector x){
    return 4*x[0]*std::pow(x[1]+36*x[0]*x[2],2)/std::pow((1+mL(x))*(1+mT(x)),2);
}

double z_der(RealVector x, double e, double d){
    return -x[3]*eta(x)+(e+d)*0;
}

double zero(double, double){
    return 0;
}


std::vector<std::function<double(RealVector, double, double)> > functions{k_der, u_der, l_der, z_der};
std::vector<std::function<double(double, double)> > scale{k_scale,u_scale,l_scale,zero};
EquationSet equations = EquationSet(functions, scale);

