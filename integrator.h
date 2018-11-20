#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <vector>
#include "realvector.h"
#include "equationset.h"


struct Plot
{
    std::vector<double> vals;
    std::vector<double> derivatives;
    std::vector<double> times;
};

std::ostream& operator << (std::ostream& out, Plot p);

class Integrator
{
public:
    Integrator();
    Integrator(EquationSet equations);

    std::vector<Plot> integrate(double start_t, double end_t, double delta_t,
                      RealVector starting_point);
    EquationSet equations;

private:
    RealVector rk4(RealVector point, double delta_t, double t, std::vector<Plot> *plots);

};

#endif // INTEGRATOR_H
