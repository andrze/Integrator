#include "integrator.h"
#include <stdexcept>
#include <cmath>
#include <ostream>

std::ostream& operator << (std::ostream& out, PlotSet p){
    for(size_t i=0; i<p.vals.size(); i++){
        out << p.times[i] << " " << p.vals[i] << std::endl;
    }

    return out;
}

Integrator::Integrator()
{

}

Integrator::Integrator(EquationSet equations){
    this->equations = equations;
}

PlotSet Integrator::integrate(double start_t, double end_t, double delta_t,
                              RealVector starting_point){

    PlotSet plots;
    double t = start_t;

    RealVector point = starting_point;
    int num_of_steps = int(std::ceil( (end_t - start_t) / delta_t));

    if(num_of_steps <= 0){
        throw std::invalid_argument( "It takes non-positive number of steps from start to end" );
    }

    plots.vals.push_back(point);
    plots.times.push_back(t);

    for(int i=0; i<num_of_steps; i++){
        point = this->rk4(point, delta_t);
        t += delta_t;

        plots.vals.push_back(point);
        plots.times.push_back(t);
    }

    return plots;
}

RealVector Integrator::rk4(RealVector point, double delta_t){

    RealVector first, second, third, fourth, final;

    first = point;

    second = first + delta_t/2*equations.evaluate(first);

    third = first + delta_t/2*equations.evaluate(second);

    fourth = first + delta_t*equations.evaluate(third);

    final = point + (equations.evaluate(first)
                     + 2*equations.evaluate(second)
                     + 2*equations.evaluate(third)
                     + equations.evaluate(fourth))*delta_t/6.;

    return final;
}




