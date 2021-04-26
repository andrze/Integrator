#include <cmath>
#include <vector>
#include <valarray>
#include <array>
#include <iostream>
#include "realvector.h"
#include "physics.h"
#include "simpson.h"
#include "regulator.h"
#include "gaussquadrature.h"

double Power(double x, double y){
    return std::pow(x,y);
}

double K_der(RealVector x, double eta, double d, double n){
    double K = x[0];
    double u20 = x[1];

    auto integrand = [&](double y){
        double y2=y*y;
        double reg = (2 - eta)*r(y2) - 2*y2*rp(y2);
        double gp = G(0, 1, y2);
        double gs = G(2*K*u20, 1, y2);

        double val = 6*Power(gs,2) + 2*Power(gp,2)*(-1 + n);

        return std::pow(y, d-1)*val*reg;
    };

    return -gauss_legendre_integrate(integrand)-(2-d-eta)*K;
}

double Z_der(RealVector x, double eta){
    return eta*x[2];
}

double find_eta(RealVector x, double d, double n){
    double K = x[0];
    double U2 = x[1];

    auto integrand = [&](double y){
        double y2 = y*y;
        double gp = G(0, 1, y2);
        double gs = G(2*K*U2, 1, y2);

        //Ising: double value = (16*Power(gp,2)*Power(gs,2)*K*Power(U2,2)*((1 + rp(y2))*(-d + 2*(gp + gs)*y2 + 2*(gp + gs)*y2*rp(y2)) - 2*y2*rp2(y2)))/d;
        double value =
                (8*K*Power(U2,2)*(36*Power(gs,5)*y2*Power(1 + rp(y2),2)
                                  + 4*Power(gp,5)*(-1 + n)*y2*Power(1 + rp(y2),2)
                                  - 9*Power(gs,4)*(d + d*rp(y2) + 2*y2*rp2(y2))
                                  - Power(gp,4)*(-1 + n)*(d + d*rp(y2) + 2*y2*rp2(y2))))/d;

        return std::pow(y, d-1)*value;
    };


    auto Z_der_num = [&](double y){

        double y2 = y*y;
        return (2*r(y2)-2*y2*rp(y2))*integrand(y);
    };

    auto Z_der_den = [&](double y){
        double y2 = y*y;

        return -r(y2)*integrand(y);
    };

    double num_integrated = gauss_legendre_integrate(Z_der_num);
    double den_integrated = gauss_legendre_integrate(Z_der_den);

    //std::cout<< num_integrated <<' '<<den_integrated<<' '<<-num_integrated/(1-den_integrated)<<'\n';
    return -num_integrated/(1-den_integrated);
}

double u20_der(RealVector x, double eta, double d, double n){
    double K = x[0];
    double U2 = x[1];

    auto integrand = [&](double y){
        double y2 = y*y;
        double reg = -(-2 + eta)*r(y2) - 2*y2*rp(y2);
        double gp = G(0, 1, y2);
        double gs = G(2*K*U2, 1, y2);

        double val = 4*(9*Power(gs,3) + Power(gp,3)*(-1 + n))*Power(U2,2);

        return std::pow(y, d-1)*val*reg;
    };

    return -gauss_legendre_integrate(integrand) - (d-4+2*eta)*U2;
}


double y_der(RealVector x, double eta, double d, double n){
    double K = x[0];
    double U2 = x[1];

    auto integrand = [&](double y){
        double y2 = y*y;
        double reg = -(-2 + eta)*r(y2) - 2*y2*rp(y2);
        double gp = G(0, 1, y2);
        double gs = G(2*K*U2, 1, y2);

        double val = 4*(9*Power(gs,3) + Power(gp,3)*(-1 + n))*Power(U2,2);

        return std::pow(y, d-1)*val*reg;
    };

    return -gauss_legendre_integrate(integrand) - (d-4+2*eta)*U2;
}


