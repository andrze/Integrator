#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <vector>
#include "realvector.h"
#include "equationset.h"
#include "plotset.h"

class Integrator
{
public:
    Integrator();
    Integrator(EquationSet equations);

    PlotSet integrate(double start_t, double end_t, double delta_t,
                      RealVector starting_point, std::vector<PlotSet >* scan=nullptr);
    EquationSet equations;

private:
    RealVector rk4(double t, PlotSet* plots);
    double d;

};

#endif // INTEGRATOR_H
