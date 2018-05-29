#include <cmath>
#include <vector>
#include "realvector.h"

const double pi = std::acos(-1);
const double d = 3., gamma3 = gamma(3);

double v_d(double d){
    return 1/(std::pow(2,d+1) * std::pow(pi,d/2) * tgamma(d/2));
}

double l(double d, int n, double x, double eta=0){
    if(n==0){
        n=1;
    }
    return (2/d)*n*(1-eta/(d+2))/(std::pow(1+x,n+1));
}

double k_der(RealVector x){
    return -1*x[0] + 6 * v_d(3.) * l(3, 1, 2*x[0]*x[1]);
}

double l_der(RealVector x){
    return -x[1] + 18 * v_d(3.) * std::pow(x[1],2) * l(3, 2, 2*x[0]*x[1]);
}

double zero(RealVector){
    return 0;
}

std::vector<std::function<double(RealVector)> > functions{k_der, l_der, zero,
                                                          zero, zero, zero};
std::vector<bool> differential{1,1,0,1,1,1};
extern const EquationSet equations(functions, differential);


