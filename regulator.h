#ifndef REGULATOR_H
#define REGULATOR_H

double yr(double y);// y r(y))
double y2r1(double y);// y^2 r'(y)
double y3r2(double y); // y^3 r''(y)
double dry(double y); //Derivative of (y r(y))
double dr1y2(double y); //Derivative of (y^2 r'(y))

double G(double m, double Z, double y);

#endif // REGULATOR_H
