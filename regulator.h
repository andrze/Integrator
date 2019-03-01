#ifndef REGULATOR_H
#define REGULATOR_H

double R(double y, double L);
double R10(double y, double L);
double R20(double y, double L);
double R01(double y, double L);
double R11(double y, double L);
double R21(double y, double L);

double G(double m, double Z, double Zp, double L, double y);

#endif // REGULATOR_H
