#include <cmath>
#include "regulator.h"

// Regulator function and derivatives divided by Z (variable regulating fluctuactions)

const double long_cutoff = 20.;
const double short_cutoff = 5e-2;

double yr(double y){
    if(y<short_cutoff){
        double y2 = y*y;
        return 2-y+y2/6-y2*y2/360;
    }
    if(y>long_cutoff){
        return 0.;
    }
    return 2*y/expm1(y);
}

double y2r1(double y){
    double y2 = y*y;
    if(y<short_cutoff){
        return -2+y2/6-y2*y2/120;
    }
    if(y>long_cutoff){
        return 0.;
    }
    double expm = expm1(y);
    double denominator = std::pow(expm, -2.);

    return -2*y2*(expm+1)*denominator;
}

double y3r2(double y){
    if(y<short_cutoff){
        double y2 = y*y;
        return 4-y2*y2/60;
    }
    if(y>long_cutoff){
        return 0.;
    }
    double expm = expm1(y);
    double expy = expm+1;
    double denominator = std::pow(expm, -3.);

    return 2*expy*(expy+1)*std::pow(y,3)*denominator;
}

double dry(double y){
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

double dr1y2(double y){

    if(y<short_cutoff){
        return y/3-std::pow(y,3)/30;
    }
    if(y>long_cutoff){
        return 0.;
    }
    double expm = expm1(y);
    double expy = expm+1;
    double denominator = std::pow(expm, -3.);

    return 2*expy*y*(-2*expm + y*(expy+1))*denominator;
}


double G(double m, double Z, double y){

    return 1/(m + Z*y + yr(y));
}
