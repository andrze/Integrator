#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <vector>
#include "realvector.h"


struct PlotSet
{
    std::vector<RealVector > vals;
    std::vector<double> times;

};

std::ostream& operator << (std::ostream& out, PlotSet p);

class Integrator
{
public:
    Integrator();
    Integrator(EquationSet equations);

    PlotSet integrate(double start_t, double end_t, double delta_t,
                      RealVector starting_point);

private:
    EquationSet equations;
    RealVector rk4(RealVector point, double delta_t);

};

#endif // INTEGRATOR_H
