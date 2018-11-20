#include "integrator.h"
#include <stdexcept>
#include <cmath>
#include <ostream>
#include <limits>
#include <iomanip>
#include <iostream>
#include "physics.h"


Integrator::Integrator()
{

}

Integrator::Integrator(EquationSet equations){
    this->equations = equations;
}

void push_point(std::vector<Plot>* plots, RealVector point, double t){
    for(size_t i=0; i < plots->size(); i++){
        (*plots)[i].vals.push_back(point[i]);
        (*plots)[i].times.push_back(-log(t));
    }
}

void push_ders(std::vector<Plot>* plots, RealVector ders){
    for(size_t i=0; i < plots->size(); i++){
        (*plots)[i].derivatives.push_back(ders[i]);
    }
}

std::vector<Plot > Integrator::integrate(double start_t, double end_t, double delta_t,
                                         RealVector starting_point){

    std::vector<Plot > plots;
    for(size_t i=0; i<starting_point.coords.size(); i++){
        plots.push_back(Plot());
    }
    double old_t = start_t;
    double t = start_t;
    double jump = exp(delta_t);

    RealVector point = starting_point;
    //int num_of_steps = int(std::ceil( (end_t - start_t) / delta_t));

    if((end_t - start_t) / delta_t <= 0){
        throw std::invalid_argument( "It takes non-positive number of steps from start to end" );
    }

    push_point(&plots, point, t);

    int n_steps=0;
    while(t>end_t){
        if(n_steps%100 == 0){
            std::cout<<-log(t)<<std::endl;
        }
        t = old_t*jump;
        double h = t-old_t;
        old_t = t;

        point = this->rk4(point, h, t, &plots);

        for(size_t j=0; j<point.coords.size(); j++){
            if(point[j]<0.){
                if(j!=2){
                    return plots;
                } else {
                    point[j] = 0.;
                }
            }
        }

        push_point(&plots, point, t);
        n_steps++;
    }

    push_ders(&plots, equations.evaluate(point, t));
    return plots;
}

RealVector Integrator::rk4(RealVector point, double delta_t, double t, std::vector<Plot> *plots){

    RealVector first=point, second, third, fourth, final;
    RealVector first_eval, second_eval, third_eval;

    first_eval = equations.evaluate(point, t);
    push_ders(plots, first_eval);

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




