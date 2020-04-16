#include <cmath>
#include "regulator.h"

// Regulator function and derivatives divided by Z (variable regulating fluctuactions)

const double long_cutoff = 20.;
const double short_cutoff = 5e-2;

double r(double y){
    if(y<short_cutoff){
        double y2 = y*y;
        return 2-y+y2/6-y2*y2/360;
    }
    if(y>long_cutoff){
        return 0.;
    }
    return 2*y/expm1(y);
}

double rp(double y){
    if(y<short_cutoff){
        return -1+y/3-std::pow(y,3)/90;
    }
    if(y>long_cutoff){
        return 0.;
    }
    double expm = expm1(y);
    double expy = expm+1;
    double denominator = std::pow(expm, -2.);

    return 2*(expm-y*expy)*denominator;
}

double rp2(double y){
    if(y<short_cutoff){
        return 1/3-y*y/30+std::pow(y,4)/504;
    }
    if(y>long_cutoff){
        return 0.;
    }
    double expm = expm1(y);
    double expy = expm+1;
    double denominator = std::pow(expm, -3.);

    return 2*expy*(-2*expm +y*(expy+1))*denominator;

}


double G(double m, double Z, double y){

    return 1/(m + Z*y + r(y));
}
