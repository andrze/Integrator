#include <cmath>
#include <vector>
#include <valarray>
#include <array>
#include <iostream>
#include "realvector.h"
#include "physics_hex.h"
#include "simpson.h"
#include "regulator.h"

double A_der_hex(RealVector x, double L, double zp_der, double d, double d_factor){
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2];      //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;
    auto integrand = [&](double y){
        return std::pow(y,d/2-1)*(x[4]*R01(L,y)+zp_der*R(L,y))*(
                    (3*u+2*Y*y*L2)*std::pow(G(x[1],x[3],x[4],L,y),2)
                     +(u+4*l)*std::pow(G(x[2],x[4],x[4],L,y),2)); };

    return x[5]*d_factor*std::pow(L,d)/(2*x[0]*x[1])*integral(integrand);
}

double Ms_der_hex(RealVector x, double L, double zp_der, double d, double d_factor){
    double a2 = x[0]*x[0];
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2]; //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;

    auto integrand = [&](double y)-> double{
        double potential = (3*u+2*Y*y*L2);
        auto Gs = G(x[1],x[3],x[4],L,y);
        auto Gp = G(x[2],x[4],x[4],L,y);

        return std::pow(y,d/2-1)*(x[4]*R01(L,y)+zp_der*R(L,y))*
                (potential*std::pow(Gs,2)+u*std::pow(Gp,2)
                 +std::pow(potential,2)*std::pow(Gs,3)
                 +std::pow(u+4*l,2)*std::pow(Gp,3));
    };

    return x[5]*d_factor*std::pow(L,d)*integral(integrand)/a2;
}

double Mp_der_hex(RealVector x, double L, double zp_der, double d, double d_factor){

    double a2 = x[0]*x[0];
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2];  //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;

    auto integrand = [&](double y)-> double{
        double uq = u+Y*L2*y;
        auto Gs = G(x[1],x[3],x[4],L,y);
        auto Gp = G(x[2],x[4],x[4],L,y);
        double Gs2 = std::pow(Gs,2), Gp2 = std::pow(Gp,2);

        return std::pow(y,d/2-1)*2*(x[4]*R01(L,y)+zp_der*R(L,y))
                    *((1/x[1])*((2*uq+3*u)*Gs2+(u+4*l)*Gp2)
                      +(11*uq+14*l)*Gs*Gp*(Gs+Gp));
    };

    return x[5]*d_factor*std::pow(L,d)*x[2]*integral(integrand)/a2;
}

double Zs_der_hex(RealVector x, double L, double zp_der, double d, double d_factor){
    double a2 = x[0]*x[0];
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2];  //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;

    auto integrand = [&](double y)-> double{
        double q2 = y*L2;
        double uq = u+Y*q2;
        double ps = u+2*uq;
        double ps2 = std::pow(ps,2);
        double pp = u+4*l;
        double pp2 = std::pow(pp,2);

        double r01 = x[4]*R01(L,y) + zp_der*R(L,y);
        double r11 = x[4]*R11(L,y) + zp_der*R10(L,y);
        double r21 = (x[4]*R21(L,y) + zp_der*R20(L,y))*q2 + r11;
        double r10 = x[4]*(1+R10(L,y));
        double sr10 = r10+Y;
        double r20 = x[4]*R20(L,y)*q2 + r10;
        double sr20 = r20+Y;

        std::array<double,6> Gs, Gp;
        Gs[1] = G(x[1],x[3],x[4],L,y);
        Gp[1] = G(x[2],x[4],x[4],L,y);
        //G[2] = ((2*uq+3*l)*std::pow(Gs,2)+ (u+6l)*std::pow(Gp,2));
        for(size_t i=2; i<6; i++){
            Gs[i] = Gs[1]*Gs[i-1];
            Gp[i] = Gp[1]*Gp[i-1];
        }

        return std::pow(y,d/2-1)*
        	   (Y/u*Gs[2]*2*uq*r01
               +Gs[3]*(0.5*ps2*r21
                      +2*r01*Y*(u+5*uq)
                      +2*q2*Y*ps*r11)
               -Gs[4]*(2*q2*ps2*sr10*r11
                       +1.5*ps2*sr20*r01
                       +6*q2*ps*sr10*r01)
               +Gs[5]*4*q2*ps2*sr10*sr10*r01
               +Y/u*Gp[2]*pp*r01
               +Gp[3]*(0.5*pp2*r21
                      +2*Y*pp*r01)
               -Gp[4]*pp2*(2*q2*r10*r11
                          +1.5*r20*r01)
               +Gp[5]*4*q2*pp2*r10*r10*r01);
    };

    return x[5]*d_factor*std::pow(L,d)*integral(integrand)/a2;
}

double Zp_der_hex(RealVector x, double L, double d, double d_factor){
    double a2 = x[0]*x[0];
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2]/3;  //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;

    auto Z_der = [&](double y){
        double q2 = y*L2;
        double potential = u+4*l+Y*y*L2;
        double p2 = std::pow(potential,2);

        double r00 = R(L,y);
        double r10 = R10(L,y);
        double r20 = R20(L,y);

        double zr10 = x[4]*(1+R10(L,y));
        double zr20 = x[4]*R20(L,y)*q2 + zr10;

        double Gs = G(x[1],x[3],x[4],L,y);
        std::array<double,5> Gp;
        Gp[0] = 1;
        Gp[1] = G(x[2],x[4],x[4],L,y);
        for(size_t i=2; i<5; i++){
            Gp[i] = std::pow(Gp[1],i);
        }

        double Z_der = (-Y*r00*Gp[2]
                        +p2*Gs*(q2*zr10*(r00*zr10*(6*Gp[4]
                                                   +2*Gp[3]*Gs)
                                         -r10*4*Gp[3])
                                +Gp[2]*(r00*zr20*(-Gs-2*Gp[1])
                                        +(q2*r20+r10))));

        return std::pow(y,d/2-1)*Z_der;
    };

    auto integrand = [&](double y){
        double q2 = y*L2;
        double potential = u+2*l+Y*y*L2;
        double p2 = std::pow(potential,2);

        double zr01 = x[4]*R01(L,y);
        double zr11 = x[4]*R11(L,y);
        double zr21 = (x[4]*R21(L,y))*q2 + zr11;
        double zr10 = x[4]*(1+R10(L,y));
        double zr20 = x[4]*R20(L,y)*q2 + zr10;

        double Gs = G(x[1],x[3],x[4],L,y);
        std::array<double,5> Gp;
        Gp[0] = 1;
        Gp[1] = G(x[2],x[4],x[4],L,y);
        for(size_t i=2; i<5; i++){
            Gp[i] = std::pow(Gp[1],i);
        }

        double Z_part = (
                        -Y*zr01*Gp[2]
                        +p2*Gs*(zr10*q2*(zr10*zr01*(6*Gp[4]+2*Gp[3]*Gs)
                                        -4*zr11*Gp[3])
                                -zr01*zr20*(Gp[2]*Gs + 2*Gp[3])
                                +zr21*Gp[2])
                            );

        return std::pow(y,d/2-1)*Z_part;
    };
    double Z_der_integrated = x[5]*d_factor*std::pow(L,d)*integral(Z_der)/a2;
    double integrated = x[5]*d_factor*std::pow(L,d)*integral(integrand)/a2;
    //std::cout<<L<<' '<<Z_der_integrated<<' '<<integrated<<' '<<integrated/(1-Z_der_integrated)<<std::endl;
    return integrated/(1-Z_der_integrated);
}

