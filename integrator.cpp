#include "integrator.h"
#include <stdexcept>
#include <cmath>
#include <ostream>
#include <limits>
#include <iomanip>
#include <iostream>


int phase_diagnosis(std::vector<Plot> result){
    if(result[0].vals.back() > 1e-03){
        double eta = -exp(-result[4].times.back())*result[4].derivatives.back()/result[4].vals.back();
        if(std::abs(eta) < 1e-04){
            return 2;
        }
    }

    for(size_t i=0; i<result[0].vals.size(); i++){
        double k = result[0].vals[i]*result[0].vals[i]*result[4].vals[i];
        if(k > 0.035 && k < 0.05){
            double eta = -exp(-result[4].times[i])*result[4].derivatives[i]/result[4].vals[i];

            if(eta > 0.26){
                return 1;
            }
        }
    }

    return 0;
}


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

void pop_point(std::vector<Plot>& plots){
    for(auto it=plots.begin(); it != plots.end(); it++){
        it->vals.pop_back();
        it->times.pop_back();
        it->derivatives.pop_back();
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
    //int contractions = 3;
    int n_steps=0;
    while(t>end_t){
        if(n_steps%500 == 0){
            std::cout<<-log(t)<<std::endl;
        }

        t = old_t*jump;
        double h = t-old_t;

        RealVector new_point = this->rk4(point, h, t, &plots, true);


        for(size_t j=0; j<point.coords.size(); j++){
            if(point[j]<0.){
                if(j!=2){
                    for(int i=0; i<-2/(10*delta_t); i++){
                        pop_point(plots);
                    }
                    std::cout<<"Zakończono w punkcie "<<point<<" na wartości "<<j<<".\n";
                    return plots;
                } else {
                    //point[j] = 0.;
                }
            }
        }
        /*bool break_for_loop = false;
        for(size_t j=0; j<point.coords.size(); j++){
            if(std::abs(plots[j].derivatives.back()*h) > 5e-3){
                std::cout<<'c';
                if(contractions == 0){
                    std::cout<<"out"<<plots[j].derivatives.size();
                    for(int i=0; i<-10; i++){
                        pop_point(plots);
                    }
                    std::cout<<"out"<<plots[j].derivatives.size();
                    return plots;
                }
                jump = std::sqrt(jump);

                t = old_t*jump;
                h = t-old_t;
                contractions--;

                new_point = this->rk4(point, h, t, &plots, false);
                break_for_loop = true;
            }
            if(break_for_loop){
                break;
            }
        }*/

        old_t = t;

        push_point(&plots, new_point, t);
        point = new_point;
        n_steps++;

        if(t<0.1 && n_steps % 10 == 0 && phase_diagnosis(plots) > 0){
            break;
        }
    }

    push_ders(&plots, equations.evaluate(point, t));
    return plots;
}

RealVector Integrator::rk4(RealVector point, double delta_t, double t, std::vector<Plot> *plots, bool add_ders=true){

    RealVector first=point, second, third, fourth, final;
    RealVector first_eval, second_eval, third_eval;

    first_eval = equations.evaluate(point, t);

    if(add_ders){
        push_ders(plots, first_eval);
    }

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




