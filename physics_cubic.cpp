#include <cmath>
#include <vector>
#include <valarray>
#include <array>
#include <iostream>
#include "realvector.h"
#include "physics_cubic.h"
#include "simpson.h"
#include "regulator.h"
#include "gaussquadrature.h"

double Power(double x, double y){
    return std::pow(x,y);
}

double K_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u30 = x[5];
    double u11 = x[6];

    auto integrand = [&](double y){
        double reg = (2 - eta)*r(y) - 2*y*rp(y);
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double val = Power(gp,2)*(1 + (2*(u01 + K*u11))/u20) +
                (Power(gs,2)*(3*u20 + 2*K*u30 + 2*y*Y))/u20;

        return std::pow(y, d/2-1)*val*reg;
    };

    return -d_factor*gauss_legendre_integrate(integrand)-(2-d-eta)*K;
}

double Z_der_cubic(RealVector x, double eta){
    return eta*x[1];
}

double find_eta(RealVector x, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u11 = x[6];

    auto integrand = [&](double y){
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double value =
                -2*Power(gp,2)*Y + (8*Power(gp,2)*gs*K*Y*
                      (d*(2*u01 + 2*K*u11 + u20) + (2 + d)*y*Y))/d +
                   (16*Power(gp,2)*Power(gs,3)*K*y*
                      Power(2*u01 + 2*K*u11 + u20 + y*Y,2)*
                      Power(1 + 2*K*Y + rp(y),2))/d +
                   Power(gs,2)*((16*Power(gp,3)*K*y*
                         Power(2*u01 + 2*K*u11 + u20 + y*Y,2)*
                         Power(1 + rp(y),2))/d -
                      (8*Power(gp,2)*K*(2*u01 + 2*K*u11 + u20 + y*Y)*
                         (4*y*Y*(1 + 2*K*Y) +
                           d*(1 + K*Y)*(2*u01 + 2*K*u11 + u20 + y*Y) +
                           (d*(2*u01 + 2*K*u11 + u20) + (4 + d)*y*Y)*
                            rp(y) + 2*y*(2*u01 + 2*K*u11 + u20 + y*Y)*
                            rp2(y)))/d);

        return std::pow(y, d/2-1)*value;
    };


    auto Z_der_num = [&](double y){

        return (2*r(y)-2*y*rp(y))*integrand(y);
    };

    auto Z_der_den = [&](double y){

        return -r(y)*integrand(y);
    };

    double num_integrated = gauss_legendre_integrate(Z_der_num)*d_factor;
    double den_integrated = gauss_legendre_integrate(Z_der_den)*d_factor;

    return -num_integrated/(1-den_integrated);
}

double Y_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u30 = x[5];
    double u11 = x[6];

    auto integrand = [&](double y){
        double reg = -(-2 + eta)*r(y) - 2*y*rp(y);
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double dimensionless =
                (Power(gp,2)*Y)/K + 4*Power(gp,3)*
                    (2*u01 + 2*K*u11 + u20)*Y -
                   (4*Power(gp,2)*gs*Y*
                      (d*(2*u01 + 2*K*u11 + u20) + (2 + d)*y*Y))/d +
                   (8*Power(gp,5)*Power(2*u01 + 2*K*u11 + u20,2)*y*
                      Power(1 + rp(y),2))/d +
                   (8*Power(gs,5)*y*Power(3*u20 + 2*K*u30 + 2*y*Y,2)*
                      Power(1 + 2*K*Y + rp(y),2))/d +
                   Power(gs,3)*(4*Y*(6*u20 + 4*K*u30 +
                         2*(2 + 1/d)*y*Y) -
                      (8*Power(gp,2)*y*
                         Power(2*u01 + 2*K*u11 + u20 + y*Y,2)*
                         Power(1 + 2*K*Y + rp(y),2))/d) -
                   (2*Power(gp,4)*Power(2*u01 + 2*K*u11 + u20,2)*
                      (d + d*rp(y) + 2*y*rp2(y)))/d -
                   (2*Power(gs,4)*(3*u20 + 2*K*u30 + 2*y*Y)*
                      ((d*(3*u20 + 2*K*u30) + 2*(4 + d)*y*Y)*
                         (1 + 2*K*Y + rp(y)) +
                        2*y*(3*u20 + 2*K*u30 + 2*y*Y)*rp2(y)))/d +
                   Power(gs,2)*(-(Y/K) -
                      (8*Power(gp,3)*y*
                         Power(2*u01 + 2*K*u11 + u20 + y*Y,2)*
                         Power(1 + rp(y),2))/d +
                      (4*Power(gp,2)*(2*u01 + 2*K*u11 + u20 + y*Y)*
                         (4*y*Y*(1 + 2*K*Y) +
                           d*(1 + K*Y)*(2*u01 + 2*K*u11 + u20 + y*Y) +
                           (d*(2*u01 + 2*K*u11 + u20) + (4 + d)*y*Y)*
                            rp(y) + 2*y*(2*u01 + 2*K*u11 + u20 + y*Y)*
                            rp2(y)))/d);

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*gauss_legendre_integrate(integrand) - (d-2+2*eta)*Y;
}

double u20_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u30 = x[5];
    double u11 = x[6];
    double u40 = x[7];
    double u21 = x[8];

    auto integrand = [&](double y){
        double reg = -(-2 + eta)*r(y) - 2*y*rp(y);
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double val = 2*Power(gp,3)*Power(2*u01 + 2*K*u11 + u20,2) +
                Power(gp,2)*(-2*(2*u11 + K*u21) +
                   (2*(u01 + K*u11)*u30)/u20) +
                2*Power(gs,3)*Power(3*u20 + 2*K*u30 + 2*y*Y,2) +
                Power(gs,2)*(-2*(u30 + K*u40) +
                   (2*u30*(K*u30 + y*Y))/u20);

        return std::pow(y, d/2-1)*val*reg;
    };

    return -d_factor*gauss_legendre_integrate(integrand) - (d-4+2*eta)*u20;
}

double u01_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u30 = x[5];
    double u11 = x[6];
    double u21 = x[8];
    double u02 = x[9];

    auto integrand = [&](double y){
        double reg = -(-2 + eta)*r(y) - 2*y*rp(y);
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double dimensionless = Power(gp,2)*(u11 + (2*u11*(u01 + K*u11))/u20) +
                2*gp*gs*(5*u11 + K*(3*u02 + u21)) +
                Power(gs,2)*((u11*(3*u20 + 2*K*u30 + 2*y*Y))/u20 -
                   4*Power(gp,2)*(y + K*(u01 + u20 + y*Y) + r(y))*
                    (-3*Power(u01,2) + 5*u11*y +
                      u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                         6*(u20 + y*Y)) +
                      Power(K,2)*(-4*Power(u11,2) +
                         6*u02*(u20 + y*Y)) +
                      K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                      (5*u11 + K*(3*u02 + u21))*r(y)));

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*gauss_legendre_integrate(integrand) - (d-4+2*eta)*u01;
}

double u30_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u30 = x[5];
    double u11 = x[6];
    double u40 = x[7];
    double u21 = x[8];

    auto integrand = [&](double y){
        double reg = -(-2 + eta)*r(y) - 2*y*rp(y);
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double dimensionless = -6*Power(gp,4)*Power(2*u01 + 2*K*u11 + u20,3) +
                6*Power(gp,3)*(2*u01 + 2*K*u11 + u20)*
                 (4*u11 + 2*K*u21 + u30) +
                Power(gp,2)*(-6*u21 + (2*(u01 + K*u11)*u40)/u20) +
                (2*Power(gs,2)*u40*(-2*u20 + K*u30 + y*Y))/u20 +
                6*Power(gs,3)*(5*u30 + 2*K*u40)*
                 (3*u20 + 2*K*u30 + 2*y*Y) -
                6*Power(gs,4)*Power(3*u20 + 2*K*u30 + 2*y*Y,3);

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*gauss_legendre_integrate(integrand) - (2*d-6+3*eta)*u30;
}

double u11_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u30 = x[5];
    double u11 = x[6];
    double u21 = x[8];
    double u02 = x[9];

    auto integrand = [&](double y){
        double reg = -(-2 + eta)*r(y) - 2*y*rp(y);
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double dimensionless =
                Power(gp,2)*(u21 + (2*(u01 + K*u11)*u21)/u20) +
                   (6*gp*(u02 + 2*u21))/(y + 2*K*(u20 + y*Y) + r(y)) +
                   8*Power(gp,2)*Power(gs,3)*(3*u20 + 2*K*u30 + 2*y*Y)*
                    (y + K*(u01 + u20 + y*Y) + r(y))*
                    (-3*Power(u01,2) + 5*u11*y +
                      u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                         6*(u20 + y*Y)) +
                      Power(K,2)*(-4*Power(u11,2) +
                         6*u02*(u20 + y*Y)) +
                      K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                      (5*u11 + K*(3*u02 + u21))*r(y)) +
                   Power(gs,2)*((u21*(3*u20 + 2*K*u30 + 2*y*Y))/u20 +
                      8*Power(gp,3)*(2*u01 + 2*K*u11 + u20)*
                       (y + K*(u01 + u20 + y*Y) + r(y))*
                       (-3*Power(u01,2) + 5*u11*y +
                         u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                            6*(u20 + y*Y)) +
                         Power(K,2)*(-4*Power(u11,2) +
                            6*u02*(u20 + y*Y)) +
                         K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                         (5*u11 + K*(3*u02 + u21))*r(y)) +
                      Power(gp,2)*(-4*(y + K*(u01 + u20 + y*Y) + r(y))*
                          (u01*(-9*u11 + K*u21 - 6*u30) +
                            Power(K,2)*(-6*u11*u21 + 6*u02*u30) +
                            y*(3*u02 + 6*u21 - 5*u11*Y) +
                            K*(-11*Power(u11,2) + 15*u02*u20 +
                               2*u20*u21 + u11*u30 + (12*u02 + u21)*y*Y
                               ) + 3*(u02 + 2*u21)*r(y)) -
                         4*(u01 + 2*u20 + K*(u11 + u30) + y*Y)*
                          (-3*Power(u01,2) + 5*u11*y +
                            u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                               6*(u20 + y*Y)) +
                            Power(K,2)*
                             (-4*Power(u11,2) + 6*u02*(u20 + y*Y)) +
                            K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                            (5*u11 + K*(3*u02 + u21))*r(y)) -
                         4*(5*u11 + K*(3*u02 + u21))*
                          (y*(u01 + 2*u20 + y*Y) +
                            2*Power(K,2)*
                             (u11*u20 + u01*u30 + u11*y*Y) +
                            K*(u20*(5*u01 + u20) + (u11 + u30)*y +
                               (4*u01 + u20)*y*Y) +
                            (u01 + 2*u20 + K*(u11 + u30) + y*Y)*r(y))));

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*gauss_legendre_integrate(integrand) - (2*d-6+3*eta)*u11;
}

double u40_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u30 = x[5];
    double u11 = x[6];
    double u40 = x[7];
    double u21 = x[8];

    auto integrand = [&](double y){
        double reg = -(-2 + eta)*r(y) - 2*y*rp(y);
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double dimensionless = 24*Power(gp,5)*Power(2*u01 + 2*K*u11 + u20,4) -
                36*Power(gp,4)*Power(2*u01 + 2*K*u11 + u20,2)*
                 (4*u11 + 2*K*u21 + u30) +
                Power(gp,3)*(6*Power(4*u11 + 2*K*u21 + u30,2) +
                   8*(2*u01 + 2*K*u11 + u20)*(6*u21 + u40)) -
                36*Power(gs,4)*(5*u30 + 2*K*u40)*
                 Power(3*u20 + 2*K*u30 + 2*y*Y,2) +
                24*Power(gs,5)*Power(3*u20 + 2*K*u30 + 2*y*Y,4) +
                Power(gs,3)*(6*Power(5*u30 + 2*K*u40,2) +
                   56*u40*(3*u20 + 2*K*u30 + 2*y*Y));

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*gauss_legendre_integrate(integrand) - (3*d-8+4*eta)*u40;
}

double u21_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u30 = x[5];
    double u11 = x[6];
    double u40 = x[7];
    double u21 = x[8];
    double u02 = x[9];

    auto integrand = [&](double y){
        double reg = -(-2 + eta)*r(y) - 2*y*rp(y);
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double dimensionless =
                -24*Power(gp,2)*Power(gs,4)*
                    Power(3*u20 + 2*K*u30 + 2*y*Y,2)*
                    (y + K*(u01 + u20 + y*Y) + r(y))*
                    (-3*Power(u01,2) + 5*u11*y +
                      u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                         6*(u20 + y*Y)) +
                      Power(K,2)*(-4*Power(u11,2) +
                         6*u02*(u20 + y*Y)) +
                      K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                      (5*u11 + K*(3*u02 + u21))*r(y)) -
                   (4*Power(gp,2)*(5*u11 + K*(3*u02 + u21))*
                      (3*Power(u20,2) + 2*u11*y + 3*u30*y +
                        2*u20*y*Y + u01*
                         (6*u20 + 9*K*u30 + 2*Power(K,2)*u40 + 4*y*Y)\
                         + 2*Power(K,2)*
                         (u20*u21 + 2*u11*u30 + u21*y*Y) +
                        K*(3*u20*u30 + y*(u21 + u40 + u30*Y) +
                           2*u11*(5*u20 + 4*y*Y)) +
                        (2*u11 + 3*u30 + K*(u21 + u40))*r(y)))/
                    Power(y + 2*K*(u20 + y*Y) + r(y),2) +
                   (16*Power(gp,3)*(5*u11 + K*(3*u02 + u21))*
                      Power(y*(u01 + 2*u20 + y*Y) +
                        2*Power(K,2)*(u11*u20 + u01*u30 + u11*y*Y) +
                        K*(u20*(5*u01 + u20) + (u11 + u30)*y +
                           (4*u01 + u20)*y*Y) +
                        (u01 + 2*u20 + K*(u11 + u30) + y*Y)*r(y),2))/
                    Power(y + 2*K*(u20 + y*Y) + r(y),3) +
                   Power(gs,3)*(-32*Power(gp,3)*
                       (2*u01 + 2*K*u11 + u20)*
                       (3*u20 + 2*K*u30 + 2*y*Y)*
                       (y + K*(u01 + u20 + y*Y) + r(y))*
                       (-3*Power(u01,2) + 5*u11*y +
                         u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                            6*(u20 + y*Y)) +
                         Power(K,2)*(-4*Power(u11,2) +
                            6*u02*(u20 + y*Y)) +
                         K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                         (5*u11 + K*(3*u02 + u21))*r(y)) +
                      4*Power(gp,2)*(4*(3*u20 + 2*K*u30 + 2*y*Y)*
                          (y + K*(u01 + u20 + y*Y) + r(y))*
                          (u01*(-9*u11 + K*u21 - 6*u30) +
                            Power(K,2)*(-6*u11*u21 + 6*u02*u30) +
                            y*(3*u02 + 6*u21 - 5*u11*Y) +
                            K*(-11*Power(u11,2) + 15*u02*u20 +
                               2*u20*u21 + u11*u30 + (12*u02 + u21)*y*Y
                               ) + 3*(u02 + 2*u21)*r(y)) +
                         4*(u01 + 2*u20 + K*(u11 + u30) + y*Y)*
                          (3*u20 + 2*K*u30 + 2*y*Y)*
                          (-3*Power(u01,2) + 5*u11*y +
                            u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                               6*(u20 + y*Y)) +
                            Power(K,2)*
                             (-4*Power(u11,2) + 6*u02*(u20 + y*Y)) +
                            K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                            (5*u11 + K*(3*u02 + u21))*r(y)) +
                         2*(5*u30 + 2*K*u40)*
                          (y + K*(u01 + u20 + y*Y) + r(y))*
                          (-3*Power(u01,2) + 5*u11*y +
                            u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                               6*(u20 + y*Y)) +
                            Power(K,2)*
                             (-4*Power(u11,2) + 6*u02*(u20 + y*Y)) +
                            K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                            (5*u11 + K*(3*u02 + u21))*r(y)))) +
                   Power(gs,2)*(-24*Power(gp,4)*
                       Power(2*u01 + 2*K*u11 + u20,2)*
                       (y + K*(u01 + u20 + y*Y) + r(y))*
                       (-3*Power(u01,2) + 5*u11*y +
                         u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                            6*(u20 + y*Y)) +
                         Power(K,2)*(-4*Power(u11,2) +
                            6*u02*(u20 + y*Y)) +
                         K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                         (5*u11 + K*(3*u02 + u21))*r(y)) +
                      Power(gp,3)*(16*(2*u01 + 2*K*u11 + u20)*
                          (y + K*(u01 + u20 + y*Y) + r(y))*
                          (u01*(-9*u11 + K*u21 - 6*u30) +
                            Power(K,2)*(-6*u11*u21 + 6*u02*u30) +
                            y*(3*u02 + 6*u21 - 5*u11*Y) +
                            K*(-11*Power(u11,2) + 15*u02*u20 +
                               2*u20*u21 + u11*u30 + (12*u02 + u21)*y*Y
                               ) + 3*(u02 + 2*u21)*r(y)) -
                         8*(4*u11 + 2*K*u21 + u30)*
                          (y + K*(u01 + u20 + y*Y) + r(y))*
                          (3*Power(u01,2) - 5*u11*y +
                            u01*(3*K*u11 - 2*Power(K,2)*u21 +
                               6*(u20 + y*Y)) +
                            Power(K,2)*
                             (4*Power(u11,2) - 6*u02*(u20 + y*Y)) -
                            K*((3*u02 + u21)*y + u11*(u20 + y*Y)) -
                            (5*u11 + K*(3*u02 + u21))*r(y)) +
                         16*(2*u01 + 2*K*u11 + u20)*
                          (u01 + 2*u20 + K*(u11 + u30) + y*Y)*
                          (-3*Power(u01,2) + 5*u11*y +
                            u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                               6*(u20 + y*Y)) +
                            Power(K,2)*
                             (-4*Power(u11,2) + 6*u02*(u20 + y*Y)) +
                            K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                            (5*u11 + K*(3*u02 + u21))*r(y))) +
                      4*Power(gp,2)*((20*Power(u11,2) +
                            u01*(8*u21 + 6*u40) +
                            u11*(33*K*u21 + 5*u30 - K*u40) +
                            u21*(-8*u20 + 6*Power(K,2)*u21 - 3*K*u30 +
                               4*y*Y) -
                            3*u02*(6*u20 + 9*K*u30 +
                               2*Power(K,2)*u40 + 4*y*Y))*
                          (y + K*(u01 + u20 + y*Y) + r(y)) -
                         2*(u01 + 2*u20 + K*(u11 + u30) + y*Y)*
                          (u01*(-9*u11 + K*u21 - 6*u30) +
                            Power(K,2)*(-6*u11*u21 + 6*u02*u30) +
                            y*(3*u02 + 6*u21 - 5*u11*Y) +
                            K*(-11*Power(u11,2) + 15*u02*u20 +
                               2*u20*u21 + u11*u30 + (12*u02 + u21)*y*Y
                               ) + 3*(u02 + 2*u21)*r(y)) +
                         (-2*u11 - 3*u30 - K*(u21 + u40))*
                          (-3*Power(u01,2) + 5*u11*y +
                            u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                               6*(u20 + y*Y)) +
                            Power(K,2)*
                             (-4*Power(u11,2) + 6*u02*(u20 + y*Y)) +
                            K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                            (5*u11 + K*(3*u02 + u21))*r(y)) -
                         6*(u02 + 2*u21)*
                          (y*(u01 + 2*u20 + y*Y) +
                            2*Power(K,2)*
                             (u11*u20 + u01*u30 + u11*y*Y) +
                            K*(u20*(5*u01 + u20) + (u11 + u30)*y +
                               (4*u01 + u20)*y*Y) +
                            (u01 + 2*u20 + K*(u11 + u30) + y*Y)*r(y))));

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*gauss_legendre_integrate(integrand) - (3*d-8+4*eta)*u21;
}

double u02_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u20 = x[3];
    double u01 = x[4];
    double u11 = x[6];
    double u21 = x[8];
    double u02 = x[9];

    auto integrand = [&](double y){
        double reg = -(-2 + eta)*r(y) - 2*y*rp(y);
        double gp = G(2*K*u01, 1, y);
        double gs = G(2*K*u20, 1+2*Y*K, y);

        double dimensionless =
                (4*Power(gp,2)*gs*(-25*Power(u11,2) +
                        2*K*u11*(u02 - u21) - 12*Power(K,2)*u02*u21 +
                        4*u01*(5*u02 + 3*u21) + 28*u02*(u20 + y*Y))*
                      (y + K*(u01 + u20 + y*Y) + r(y)))/
                    (y + 2*K*(u20 + y*Y) + r(y)) +
                   Power(gs,2)*(8*Power(gp,2)*
                       (-5*u11 - K*(3*u02 + u21))*
                       (-3*Power(u01,2) + 5*u11*y +
                         u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                            6*(u20 + y*Y)) +
                         Power(K,2)*(-4*Power(u11,2) +
                            6*u02*(u20 + y*Y)) +
                         K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                         (5*u11 + K*(3*u02 + u21))*r(y)) +
                      (16*Power(gp,3)*(y + K*(u01 + u20 + y*Y) + r(y))*
                         Power(-3*Power(u01,2) + 5*u11*y +
                           u01*(-3*K*u11 + 2*Power(K,2)*u21 -
                              6*(u20 + y*Y)) +
                           Power(K,2)*(-4*Power(u11,2) +
                              6*u02*(u20 + y*Y)) +
                           K*((3*u02 + u21)*y + u11*(u20 + y*Y)) +
                           (5*u11 + K*(3*u02 + u21))*r(y),2))/
                       (y + 2*K*(u20 + y*Y) + r(y)));

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*gauss_legendre_integrate(integrand) - (3*d-8+4*eta)*u02;
}

