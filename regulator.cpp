#include <cmath>
#include "regulator.h"

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

double R10(double, double y){
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
