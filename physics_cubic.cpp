#include <cmath>
#include <vector>
#include <valarray>
#include <array>
#include <iostream>
#include "realvector.h"
#include "physics_cubic.h"
#include "simpson.h"
#include "regulator.h"


double A_der_cubic(RealVector x, double L, double zp_der, double d, double d_factor){
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2];      //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;
    auto integrand = [&](double Q){
        double y = std::pow(Q, 2/d);
        double dimensionless = (x[4]*R01(y,L)+zp_der*R(y,L))*(
                    (3*u+2*Y*y*L2)*std::pow(G(x[1],x[3],x[4],L,y),2)
                     +(u+2*l)*std::pow(G(x[2],x[4],x[4],L,y),2));

        return dimensionless;
    };

    return x[5]*d_factor*std::pow(L,d)*integral(integrand)/(2*x[0]*x[1]);
}

double Ms_der_cubic(RealVector x, double L, double zp_der, double d, double d_factor){
    double a2 = x[0]*x[0];
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2]; //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;

    auto integrand = [&](double Q){
        double y = std::pow(Q, 2/d);

        double s_vertex = (3*u+2*Y*y*L2);
        double p_vertex = u+2*l;
        auto Gs = G(x[1],x[3],x[4],L,y);
        auto Gp = G(x[2],x[4],x[4],L,y);

        double dimensionless = (x[4]*R01(y,L)+zp_der*R(y,L))*
                (s_vertex*std::pow(Gs,2)+p_vertex*std::pow(Gp,2)
                 +std::pow(s_vertex,2)*std::pow(Gs,3)
                 +std::pow(p_vertex,2)*std::pow(Gp,3));

        return dimensionless;
    };

    return x[5]*d_factor*std::pow(L,d)*integral(integrand)/a2;
}

double Mp_der_cubic(RealVector x, double L, double zp_der, double d, double d_factor){

    double a2 = x[0]*x[0];
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2];  //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;

    auto integrand = [&](double Q)-> double{
        double y = std::pow(Q, 2/d);
        double uq = u+Y*L2*y;
        auto Gs = G(x[1],x[3],x[4],L,y);
        auto Gp = G(x[2],x[4],x[4],L,y);
        double Gs2 = std::pow(Gs,2);
        double Gp2 = std::pow(Gp,2);


        double dimensionless = (x[4]*R01(y,L)+zp_der*R(y,L))
                   *((1/x[1])*((u+2*uq)*Gs2+(u+2*l)*Gp2)
                     +(6*uq+3*l)*Gs*Gp*(Gs+Gp));

        return dimensionless;
    };

    return x[5]*d_factor*std::pow(L,d)*x[2]*integral(integrand)/a2;
}

double Zs_der_cubic(RealVector x, double L, double zp_der, double d, double d_factor){
    double a2 = x[0]*x[0];
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2];  //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;
    double Z = x[4];

    auto integrand = [&](double Q)-> double{
        double y = std::pow(Q, 2/d);
        double q2 = y*L2*2/d;
        double uq = u+Y*q2;
        double s_vertex = u+2*uq;
        double s_vertex2 = std::pow(s_vertex,2);
        double p_vertex = u+2*l;
        double p_vertex2 = std::pow(p_vertex,2);

        double r01 = Z*R01(y,L) + zp_der*R(y,L);
        double r11 = Z*R11(y,L) + zp_der*R10(y,L);
        double r21 = (Z*R21(y,L) + zp_der*R20(y,L))*q2 + r11;
        double r10 = Z*(1+R10(y,L));
        double sr10 = r10+Y;
        double r20 = Z*R20(y,L)*q2 + r10;
        double sr20 = r20+Y;

        std::array<double,6> Gs, Gp;
        Gs[1] = G(x[1],x[3],x[4],L,y);
        Gp[1] = G(x[2],x[4],x[4],L,y);
        for(size_t i=2; i<6; i++){
            Gs[i] = Gs[1]*Gs[i-1];
            Gp[i] = Gp[1]*Gp[i-1];
        }

        double dimensionless =
        	   (Y/u*Gs[2]*2*uq*r01
               +Gs[3]*(0.5*s_vertex2*r21
                      +2*r01*Y*(2*u+4*uq+q2*Y)
                      +2*q2*Y*s_vertex*r11)
               -Gs[4]*(2*q2*s_vertex2*sr10*r11
                       +1.5*s_vertex2*sr20*r01
                       +6*q2*s_vertex*sr10*r01)
               +Gs[5]*4*q2*s_vertex2*sr10*sr10*r01
               +Y/u*Gp[2]*p_vertex*r01
               +Gp[3]*(0.5*p_vertex2*r21
                      +2*Y*p_vertex*r01)
               -Gp[4]*p_vertex2*(2*q2*r10*r11
                                +1.5*r20*r01)
               +Gp[5]*4*q2*p_vertex2*r10*r10*r01);

        return dimensionless;
    };

    return x[5]*d_factor*std::pow(L,d)*integral(integrand)/a2;
}

double Zp_der_cubic(RealVector x, double L, double d, double d_factor){
    double a2 = x[0]*x[0];
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2];  //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;
    double Z = x[4];

    auto Z_der = [&](double Q){
        double y = std::pow(Q, 2/d);
        double q2 = y*L2*2/d;
        double vertex = u+2*l+Y*y*L2;
        double vertex2 = std::pow(vertex,2);

        double r00 = R(y,L);
        double r10 = R10(y,L);
        double r20 = q2*R20(y,L)+r10;

        double zr10 = Z*(1+R10(y,L));
        double zr20 = Z*R20(y,L)*q2 + zr10;

        double Gs = G(x[1],x[3],x[4],L,y);
        std::array<double,5> Gp;
        Gp[0] = 1;
        Gp[1] = G(x[2],x[4],x[4],L,y);
        for(size_t i=2; i<5; i++){
            Gp[i] = std::pow(Gp[1],i);
        }

        double Z_der = -Y*r00*Gp[2]
                       +vertex2*Gs*(q2*zr10*(r00*zr10*(6*Gp[4]
                                                      +2*Gp[3]*Gs)
                                            -r10*4*Gp[3])
                                   +Gp[2]*(r00*zr20*(-Gs-2*Gp[1])+r20));

        return Z_der;
    };

    auto integrand = [&](double Q){
        double y = std::pow(Q, 2/d);
        double q2 = y*L2*2/d;
        double vertex = u+2*l+Y*y*L2;
        double vertex2 = std::pow(vertex,2);

        double zr01 = Z*R01(y,L);
        double zr11 = Z*R11(y,L);
        double zr21 = (Z*R21(y,L))*q2 + zr11;
        double zr10 = Z*(1+R10(y,L));
        double zr20 = Z*R20(y,L)*q2 + zr10;

        double Gs = G(x[1],x[3],x[4],L,y);
        std::array<double,5> Gp;
        Gp[0] = 1;
        Gp[1] = G(x[2],x[4],x[4],L,y);
        for(size_t i=2; i<5; i++){
            Gp[i] = std::pow(Gp[1],i);
        }

        double Z_part = -Y*zr01*Gp[2]
                        +vertex2*Gs*(q2*zr10*(zr01*zr10*(6*Gp[4]
                                                        +2*Gp[3]*Gs)
                                             -zr11*4*Gp[3])
                                    +Gp[2]*(zr01*zr20*(-Gs-2*Gp[1])+zr21));

        return Z_part;
    };

    double Z_der_integrated = x[5]*d_factor*std::pow(L,d)*integral(Z_der)/a2;
    double integrated = x[5]*d_factor*std::pow(L,d)*integral(integrand)/a2;

    //std::cout<<L<<' '<<Z_der_integrated<<' '<<integrated<<' '<<integrated/(1-Z_der_integrated)<<std::endl;
    return integrated/(1-Z_der_integrated);
}

