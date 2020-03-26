#include <cmath>
#include <vector>
#include <valarray>
#include <array>
#include <iostream>
#include "realvector.h"
#include "physics_cubic.h"
#include "simpson.h"
#include "regulator.h"

double K_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u = x[3];
    double l = x[4];
    double v = x[5];
    double J = x[6];

    auto integrand = [&](double y){
        double Gp = G(2*K*l, 1, y);
        double Gs = G(2*K*u, 1 + 2*K*Y, y);
        double Gp2 = Gp*Gp, Gs2 = Gs*Gs;
        double reg = (eta*yr(y) + 2*y2r1(y));

        double dimensionless =
               -((2*l+u+2*K*J)*Gp2/u + Gs2*(3 + 2*(y*Y+ K*v)/u))*reg;

        return std::pow(y, d/2-1)*dimensionless;
    };

    return -d_factor*integral(integrand)-(2-d-eta)*K;
}

double u_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u = x[3];
    double l = x[4];
    double v = x[5];
    double J = x[6];

    auto integrand = [&](double y){
        double Gp = G(2*K*l, 1, y);
        double Gs = G(2*K*u, 1 + 2*K*Y, y);
        double Gp2 = Gp*Gp, Gs2 = Gs*Gs;
        double Gp3 = Gp2*Gp, Gs3 = Gs2*Gs;
        double reg = (eta*yr(y) + 2*y2r1(y));

        double dimensionless = -2*Gp3*std::pow(2*l + 2*K*J + u,2.) - 2*Gs3*std::pow(2*y*Y + 3*u + 2*K*v,2)
                + Gp2*(4*J-2*(l+K*J)*v/u) + 2* Gs2*v*(u-K*v-y*Y)/u ;


        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*integral(integrand) - (d-4+2*eta)*u;
}

double l_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u = x[3];
    double l = x[4];
    double J = x[6];

    auto integrand = [&](double y)-> double{
        double Gp = G(2*K*l, 1, y);
        double Gs = G(2*K*u, 1 + 2*K*Y, y);
        double Gp2 = Gp*Gp, Gs2 = Gs*Gs;

        double reg = (eta*yr(y) + 2*y2r1(y));

        double dimensionless = -2*Gp2*J*(l+K*J-2*u)/u
                +2*Gs2*J*(u-K*J-y*Y)/u
                -2*Gp*Gs*(Gp+Gs)*(3*l+2*K*J)*(l+2*y*Y+2*K*J+2*u);

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*integral(integrand) - (d-4+2*eta)*l;
}

double v_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u = x[3];
    double l = x[4];
    double v = x[5];
    double J = x[6];

    auto integrand = [&](double y){
        double Gp = G(2*K*l, 1, y);
        double Gs = G(2*K*u, 1 + 2*K*Y, y);
        double Gp2 = Gp*Gp, Gs2 = Gs*Gs;
        double Gp3 = Gp2*Gp, Gs3 = Gs2*Gs;
        double Gp4 = Gp3*Gp, Gs4 = Gs3*Gs;
        double reg = (eta*yr(y) + 2*y2r1(y));
        double V3s = 2*y*Y + 3*u + 2*K*v, V3p = 2*l + u + 2*K*J;

        double dimensionless = 6*Gp4*std::pow(V3p,3.) + 6*Gs4*std::pow(V3s,3.)
                - 6*Gp3*(v+4*J)*V3p - 30* Gs3*v*V3s;

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*integral(integrand) - (2*d-6+3*eta)*v;
}


double Y_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u = x[3];
    double l = x[4];
    double v = x[5];
    double J = x[6];

    double V3p = 2*l + u + 2*K*J;

    auto integrand = [&](double y){
        double Gp = G(2*K*l, 1, y);
        double Gs = G(2*K*u, 1 + 2*K*Y, y);
        double Gp2 = Gp*Gp, Gs2 = Gs*Gs;
        double Gp3 = Gp2*Gp, Gs3 = Gs2*Gs;
        double Gp4 = Gp3*Gp, Gs4 = Gs3*Gs;
        double Gp5 = Gp4*Gp, Gs5 = Gs4*Gs;
        double reg = (eta*yr(y) + 2*y2r1(y));
        double V3s = 2*y*Y + 3*u + 2*K*v, V3pY = V3p + y*Y;

        double rp = 1 + dry(y), rs = rp + 2*Y*K;
        double rp2 = rp*rp, rs2 = rs*rs;
        double V3s2 = V3s*V3s, V3p2 = V3p*V3p;

        double dimensionless =
               d*Y/K * (Gs2-Gp2)

               + 4*Gs*Gp2*Y*(2*y*Y + d*V3pY)
               - 4*Gp3*Y*d*V3p
               - 8*Gs3*Y*(y*Y + d*V3s)

               - 4*Gs2*Gp2*V3pY
                    *(V3pY*(d*rp + 2*dr1y2(y) + d*K*Y)
                        + 4*y*Y*rs)

               + 2*Gp4*V3p2*(d*rp+2*dr1y2(y))
               + 2*Gs4*V3s*(8*y*Y*rs + d*V3s*rs + 2*V3s*dr1y2(y))

               + 8*y*Gs2*Gp2*V3pY*V3pY*(Gs*rs2+Gp*rp2)

               - 8*Gs5*y*V3s2*rs2
               - 8*Gp5*y*V3p2*rp2;


        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor/d*integral(integrand) - (d-2+2*eta)*Y;
}

double Z_der_cubic(RealVector x, double eta){
    return eta*x[1];
}

double J_der_cubic(RealVector x, double eta, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u = x[3];
    double l = x[4];
    double v = x[5];
    double J = x[6];

    auto integrand = [&](double y){
        double Gp = G(2*K*l, 1, y);
        double Gs = G(2*K*u, 1 + 2*K*Y, y);
        double Gp2 = Gp*Gp, Gs2 = Gs*Gs;
        double Gp3 = Gp2*Gp, Gs3 = Gs2*Gs;
        double reg = (eta*yr(y) + 2*y2r1(y));
        double V3s = 2*y*Y + 3*u + 2*K*v, V3p = 2*l + u + 2*K*J;

        double dimensionless =
                - 10*J*(Gp3*V3p+Gs3*V3s)
                +2*Gp2*Gs*(v*(K*J-6*l)-J*(19*l+21*K*J+5*(y*Y+u)))
                -2*Gs2*Gp*(11*K*J*J+6*l*v+3*J*(5*y*Y+3*l+5*u+3*K*v))
                +4*Gp*Gs*(Gs2*V3s+Gp2*V3p)*(3*l+2*K*J)*(l+2*(y*Y+K*J+u))
                +4*Gp2*Gs2*(3*l*l*(l+3*y*Y+2*K*J+4*u+K*v)
                           +K*J*(4*K*K*J*J-(y*Y+u)*(y*Y-3*u+K*v)+K*J*(13*y*Y+17*u+4*K*v))
                           +l*(7*K*K*J*J+6*(y*Y+u)*(y*Y+2*u+K*v)+K*J*(28*y*Y+36*u+13*K*v)));

        return std::pow(y, d/2-1)*dimensionless*reg;
    };

    return -d_factor*integral(integrand) - (2*d-6+3*eta)*J;
}

double find_eta(RealVector x, double d, double d_factor){
    double K = x[0];
    double Y = x[2];
    double u = x[3];
    double l = x[4];
    double J = x[6];

    auto integrand = [&](double y){
        double Gp = G(2*K*l, 1, y);
        double Gs = G(2*K*u, 1 + 2*K*Y, y);
        double Gp2 = Gp*Gp, Gs2 = Gs*Gs;

        double V3pY = 2*l + 2*K*J + u + y*Y;
        double rp = 1 + dry(y), rs = 2*Y*K + rp;

        double value = 2*d*Gp2*Y
                 - 8*Gp2*Gs*Y*K*(d*V3pY + 2*y*Y)
                 - 16*Gp2*Gs2*y*K*std::pow(V3pY,2)
                        *(Gs*std::pow(rs, 2) + Gp*std::pow(rp, 2))
                 + 8*Gs2*Gp2*K*V3pY*(V3pY*(d*rp + 2*dr1y2(y) + d*K*Y)
                                     + 4*y*Y*rs);

        return std::pow(y, d/2-1)*value;
    };


    auto Z_der_num = [&](double y){

        return 2*y2r1(y)*integrand(y);
    };

    auto Z_der_den = [&](double y){

        return yr(y)*integrand(y);
    };

    double num_integrated = integral(Z_der_num)*d_factor/d;
    double den_integrated = integral(Z_der_den)*d_factor/d;

    return -num_integrated/(1-den_integrated);
}

