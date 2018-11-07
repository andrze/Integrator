#include "integrator.h"
#include <stdexcept>
#include <cmath>
#include <ostream>
#include <limits>
#include <iomanip>
#include <iostream>
#include "physics.h"

std::vector<std::vector<double> > PlotSet::transpose(){

    std::vector<std::vector<double> > result;

    for(size_t j=0; j<vals[0].coords.size(); j++){
        std::vector<double> line;
        for(size_t i=0; i<vals.size(); i++){
            line.push_back(vals[i][j]);
        }
        result.push_back(line);
    }

    return result;
}

std::vector<double> PlotSet::time_exp(){
    std::vector<double> new_times;

    for(auto t: times){
        new_times.push_back(exp(t));
    }

    return new_times;
}

std::ostream& operator << (std::ostream& out, PlotSet p){
    out << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
    if(p.times.empty()){
        for(size_t i=0; i<p.vals.size(); i++){
            out << p.vals[i] << std::endl;
        }
    } else {
        for(size_t i=0; i<p.vals.size(); i++){
            out << p.times[i] << " " << p.vals[i] << std::endl;
        }
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
    //int num_of_steps = int(std::ceil( (end_t - start_t) / delta_t));

    if((end_t - start_t) / delta_t <= 0){
        throw std::invalid_argument( "It takes non-positive number of steps from start to end" );
    }

    plots.vals.push_back(point);
    plots.times.push_back(t);

    int n_steps=0;
    while(t>end_t){
        if(n_steps%100 == 0){
            std::cout<<-log(t)<<std::endl;
        }
        double h = delta_t*t/start_t;
        point = this->rk4(point, h, t, &plots);
        t += h;

        for(size_t j=0; j<point.coords.size(); j++){
            if(point[j]<0.){
                if(j!=2){
                    plots.derivatives.push_back(equations.evaluate(plots.vals.back(),
                                                                   plots.times.back()));
                    return plots;
                } else {
                    point[j] = 0.;
                }
            }
        }

        plots.vals.push_back(point);
        plots.times.push_back(t);
        n_steps++;
    }

    plots.derivatives.push_back(equations.evaluate(plots.vals.back(),
                                                   plots.times.back()));
    return plots;
}

RealVector Integrator::rk4(RealVector point, double delta_t, double t, PlotSet *plots){

    RealVector first=point, second, third, fourth, final;
    RealVector first_eval, second_eval, third_eval;

    first_eval = equations.evaluate(point, t);
    plots->derivatives.push_back(first_eval);

    second = first + delta_t/2*first_eval;
    second_eval = equations.evaluate(second, t+delta_t/2);

    third = first + delta_t/2*second_eval;
    third_eval = equations.evaluate(third, t+delta_t/2);

    fourth = first + delta_t*third_eval;

    final = first + (first_eval
                     + 2*second_eval
                     + 2*third_eval
                     + equations.evaluate(fourth, t+delta_t))*delta_t/6.;

    return final;
}




