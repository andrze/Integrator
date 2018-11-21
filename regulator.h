#ifndef REGULATOR_H
#define REGULATOR_H

double R(double L, double y);
double R10(double L, double y);
double R20(double L, double y);
double R01(double L, double y);
double R11(double L, double y);
double R21(double L, double y);

double G(double m, double Z, double Zp, double L, double y);

#endif // REGULATOR_H
