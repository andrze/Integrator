#include <cmath>
#include "regulator.h"

// Regulator function and derivatives divided by Z (variable regulating fluctuactions)

const double long_cutoff = 20.;
const double short_cutoff = 5e-2;

double R(double y, double L){
    if(y<short_cutoff){
        double y2 = y*y;
        return 2*L*L*(1-y/2+y2/12-y2*y2/720);
    }
    if(y>long_cutoff){
        return 0.;
    }
    return 2*y*L*L/expm1(y);
}

double R10(double y, double){
    if(y<short_cutoff){
        return -(1-y/3+y*y*y/90);
    }
    if(y>long_cutoff){
        return 0.;
    }
    double expm = expm1(y);
    double denominator = std::pow(expm, -2.);

    return 2*(expm - y*(expm+1))*denominator;
}

double R20(double y, double L){
    if(y<short_cutoff){
        double y2 = y*y;
        return (1-y2/10+y2*y2/168)/(3*L*L);
    }
    if(y>long_cutoff){
        return 0.;
    }
    double expm = expm1(y);
    double expy = expm+1;
    double denominator = std::pow(expm, -3.);

    return 2*expy/(L*L)*(2 + y - (2-y)*expy)*denominator;
}

double R01(double y, double L){
    double y2 = y*y;
    if(y<short_cutoff){
        return L*(4-y2/3+y2*y2/60);
    }
    if(y>long_cutoff){
        return 0.;
    }
    double expm = expm1(y);
    double denominator = std::pow(expm,-2.);

    return 4*y2*L*(expm+1)*denominator;
}

double R11(double y, double L){
    if(y<short_cutoff){
        return -1/L*(y*2/3-y*y*y/15);
    }
    if(y>long_cutoff){
        return 0.;
    }
    double expm = expm1(y);
    double denominator = std::pow(expm,-3.);
    double expy = expm+1;

    return -4*y/L*expy*((y-2)*expy+2+y)*denominator;
}

double R21(double y, double L){
    if(y<short_cutoff){
        double y2 = y*y;
        return -std::pow(L,-3.)*(2./3-y2/5+5*y2*y2/252);
    }
    if(y>long_cutoff){
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

    return 1/(m + Z*L*L*y + Zp*R(y, L));
}
