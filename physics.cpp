#include <cmath>
#include <vector>
#include <valarray>
#include <array>
#include <iostream>
#include "realvector.h"
#include "physics.h"
#include "simpson.cpp"

const double pi = std::acos(-1);

double R(double L, double y){
    if(y<5e-2){
        double y2 = y*y;
        return 2*L*L*(1-y/2+y2/12-y2*y2/720);
    }
    if(y>20){
        return 0.;
    }
    return 2*y*L*L/expm1(y);
}

double R10(double , double y){
    if(y<5e-2){
        return -(1-y/3+y*y*y/90);
    }
    if(y>20){
        return 0.;
    }
    double expm = expm1(y);
    double denominator = std::pow(expm, -2.);

    return 2*(expm - y*(expm+1))*denominator;
}

double R20(double L, double y){
    if(y<5e-2){
        double y2 = y*y;
        return (1-y2/10+y2*y2/168)/(3*L*L);
    }
    if(y>20){
        return 0.;
    }
    double expm = expm1(y);
    double expy = expm+1;
    double denominator = std::pow(expm, -3.);

    return 2*expy/(L*L)*(2 + y - (2-y)*expy)*denominator;
}

double R01(double L, double y){
    double y2 = y*y;
    if(y<5e-2){
        return 2*L*(2-y+y2/6-y2*y2/360);
    }
    if(y>20){
        return 0.;
    }
    double expm = expm1(y);
    double denominator = std::pow(expm,-2.);

    return 4*y2*L*(expm+1)*denominator;
}

double R11(double L, double y){
    if(y<5e-2){
        return -1/L*(y*2/3-y*y*y/15);
    }
    if(y>20){
        return 0.;
    }
    double expm = expm1(y);
    double denominator = std::pow(expm,-3.);
    double expy = expm+1;

    return -4*y/L*expy*((y-2)*expy+2+y)*denominator;
}

double R21(double L, double y){
    if(y<5e-2){
        double y2 = y*y;
        return -std::pow(L,-3.)*(2-y2/5+5*y2*y2/252);
    }
    if(y>20){
        return 0.;
    }
    double expm = expm1(y);
    double denominator = std::pow(expm,-4);
    double expy = expm+1;
    double exp2y = expy*expy;
    double y2 = y*y;

    return 4*expy*std::pow(L,-3.)*(2+4*y+y2-4*(1-y2)*expy+(2-4*y+y2)*exp2y)*denominator;
}

double G(double m, double Z, double Zp, double L, double y){

    return 1/(m + Z*L*L*y + Zp*R(L, y));
}

double A_der(RealVector x, double L, double zp_der){
    if(x[0] <= 1e-12){
        return 0;
    }
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2];      //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;
    auto integrand = [&](double y){
        return (x[4]*R01(L,y)+zp_der*R(L,y))*(
                    (3*u+2*Y*y*L2)*std::pow(G(x[1],x[3],x[4],L,y),2)
                     +(u+2*l)*std::pow(G(x[2],x[4],x[4],L,y),2)); };

    return x[5]*L2/(8*pi*x[0]*x[1])*integral(integrand);
}

double Ms_der(RealVector x, double L, double zp_der){
    if(x[0] <= 1e-12){
        return 0;
    }
    double a2 = x[0]*x[0];
    double u = x[1]; //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2]; //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;

    auto integrand = [&](double y)-> double{
        double potential = (3*u+2*Y*y*L2);
        auto Gs = G(x[1],x[3],x[4],L,y);
        auto Gp = G(x[2],x[4],x[4],L,y);

        return (x[4]*R01(L,y)+zp_der*R(L,y))*
                (potential*std::pow(Gs,2)+(u+2*l)*std::pow(Gp,2)
                 +std::pow(potential,2)*std::pow(Gs,3)
                 +std::pow(u+2*l,2)*std::pow(Gp,3));
    };

    return x[5]*L2/(4*pi*a2)*integral(integrand);
}

double Mp_der(RealVector x, double L, double zp_der){

    /*if(x[0] <= 1e-6){
        return 0;
    }*/
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

        return (x[4]*R01(L,y)+zp_der*R(L,y))
                   *((1/x[1])*((2*uq+u)*Gs2+(u+2*l)*Gp2)
                     +(6*uq+3*l)*Gs*Gp*(Gs+Gp));
    };

    return x[5]*L2*x[2]/(4*pi*a2)*integral(integrand);
}

double Zs_der(RealVector x, double L, double zp_der){
    if(x[0] <= 1e-6){
        return 0;
    }
    double a2 = x[0]*x[0];
    double u = (x[1]-x[2]); //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2];  //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;

    auto integrand = [&](double y)-> double{
        double q2 = y*L2;
        double uq = u+Y*q2;
        double ps = u+2*uq;
        double ps2 = std::pow(ps,2);
        double pp = u+2*l;
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

        return Y/u*Gs[2]*2*uq*r01
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
               +Gp[5]*4*q2*pp2*r10*r10*r01;
    };

    return x[5]*L2/(4*pi*a2)*integral(integrand);
}

double Zp_der(RealVector x, double L){
    if(x[0] <= 1e-6){
        return 0;
    }
    double a2 = x[0]*x[0];
    double u = (x[1]-x[2]); //Wartości u, l i Y są tu pomnożone przez alpha2
    double l = x[2]/3;  //Żeby uniknąć osobliwości
    double Y = (x[3]-x[4]);
    double L2 = L*L;

    auto Z_der = [&](double y){
        double q2 = y*L2;
        double potential = u+2*l+Y*y*L2;
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

        return Z_der;
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

        return Z_part;
    };
    double Z_der_integrated = x[5]*L2/(4*pi*a2)*integral(Z_der);
    double integrated = x[5]*L2/(4*pi*a2)*integral(integrand);
    //std::cout<<L<<' '<<Z_der_integrated<<' '<<integrated<<' '<<integrated/(1-Z_der_integrated)<<std::endl;
    return integrated/(1-Z_der_integrated);
}

