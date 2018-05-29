#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <vector>
#include "realvector.h"


struct PlotSet
{
    std::vector<RealVector > vals;
    std::vector<double> times;

    std::vector<std::vector<double> > transpose();
    std::vector<double> time_exp();
};

std::ostream& operator << (std::ostream& out, PlotSet p);

class Integrator
{
public:
    Integrator();
    Integrator(EquationSet equations);

    PlotSet integrate(double start_t, double end_t, double delta_t,
                      RealVector starting_point);
    EquationSet equations;

private:
    RealVector rk4(RealVector point, double delta_t);

};

#endif // INTEGRATOR_H
