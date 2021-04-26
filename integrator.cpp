#include "integrator.h"
#include <stdexcept>
#include <cmath>
#include <ostream>
#include <limits>
#include <iomanip>
#include <iostream>
#include <set>

int signum(double x){
    if(x>0){
        return 1;
    } else if(x<0){
        return -1;
    }
    return 0;
}

Integrator::Integrator()
{
    this->d = -1;
}

Integrator::Integrator(EquationSet equations){
    this->equations = equations;
    this->d = equations.d;
}

PlotSet Integrator::integrate(double start_t, double end_t, double delta_t,
                              RealVector starting_point, std::vector<PlotSet >* scan){

    PlotSet plots(starting_point.coords.size(), d);
    double old_t = start_t;
    double t = start_t;
    int ordered_phase_stalling = 5;
    if(scan){
        scan->push_back(PlotSet(starting_point.coords.size(), d));
    }

    RealVector point = starting_point;

    if((end_t - start_t) / delta_t <= 0){
        throw std::invalid_argument( "It takes non-positive number of steps from start to end" );
    }

    int n_steps=0;
    while(t < end_t){
        if(n_steps%500 == 0){
            std::cout<<std::setprecision(3);
            std::cout<< t <<", "<<std::flush;
            std::cout<<std::setprecision(8);
        }
        auto ders = equations.evaluate(point);
        plots.push_to_each(point, ders, t);

        RealVector new_point = this->rk4(delta_t, &plots);
        if(new_point[0] <= 0.){
            plots.phase=0;
            return plots;
        }

        bool discontinuity = false;
        for(size_t j=0; j<point.coords.size(); j++){
            bool parameter_jump = false, parameter_changed_sign = false;
            if(std::abs(point[j]) <= std::numeric_limits<double>::epsilon()) continue;

            if(j==3){
                parameter_changed_sign = (signum(new_point[j]) != signum(point[j]));
            }
            if(j==0||j==1||j==5){
                parameter_jump =  std::abs(std::log(std::abs(new_point[j]/point[j]))) > 0.5;
            }
            if(parameter_changed_sign || parameter_jump){
                discontinuity = true;
            }
        }
        if(discontinuity){
            plots.pop_from_each();
            if(delta_t<1e-3){
                auto rescaled_kappa = plots.rescaled(0, size_t(n_steps-1));
                double log_der = std::abs(rescaled_kappa.second/rescaled_kappa.first);
                if(log_der>-0.1){
                    plots.phase = 0;
                    std::cout<<"Zakończono w fazie nieuporządkowanej (nieciągłość kappa)\n";
                } else {
                    plots.phase = 2;
                    std::cout<<"Zakończono w fazie uporządkowanej (nieciągłość potencjału)\n";
                }
                return plots;
            }
            delta_t /= 2.;
            continue;
        }

        old_t = t;
        t = old_t + delta_t;

        if( static_cast<int>(t) != static_cast<int>(old_t)){
            if(scan && static_cast<int>(t) <= 50){
                scan->back().push_to_each(plots.back_vals(), plots.back_ders(), plots.back_time());
            }
        }

        point = new_point;
        n_steps++;

        if(n_steps % 3 == 0 && plots.phase_diagnosis() == 2){
            if(ordered_phase_stalling>0){
                ordered_phase_stalling--;
            } else {
                std::cout<<"Zakończono w fazie uporządkowanej\n";
                plots.phase = 2;
                return plots;
            }
        }
    }
    std::cerr<<"Zakończono na końcu czasu symulacji\n";
    plots.phase = -1;
    return plots;
}

RealVector Integrator::rk4(double delta_t, PlotSet* plots){

    RealVector first=plots->back_vals(), second, third, fourth, final;
    RealVector first_eval, second_eval, third_eval;


    first_eval = plots->back_ders();

    second = first + delta_t/2*first_eval;
    second_eval = equations.evaluate(second);

    third = first + delta_t/2*second_eval;
    third_eval = equations.evaluate(third);

    fourth = first + delta_t*third_eval;

    final = first + (first_eval
                     + 2*second_eval
                     + 2*third_eval
                     + equations.evaluate(fourth))*delta_t/6.;

    return final;
}




