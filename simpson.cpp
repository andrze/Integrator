#include <cmath>
#include <iostream>
#include <vector>
#include <functional>
#include <limits>

double integral(std::function<double(double)> func){
    double h0=0.01;
    double q=0.;
    double result=0.;
    int num_steps=51;
    double f0 = func(0.);

    for(int i=0; i<100; i++){
        double coef = std::pow(std::abs(func(q)/f0),0.2);
        double h;
        if(coef >= h0*std::numeric_limits<double>::epsilon()){
            h = h0/coef;
        } else {
            h = h0;
        }
        double partial_result=func(q);
        for(int j=0; j<num_steps; j++){
            q += h;
            double f = func(q);
            if(j % 2 == 1){
                partial_result += 2*f;
            } else {
                partial_result += 4*f;
            }
        }
        q+=h;
        partial_result += func(q);
        result += partial_result*h/3;

        if(std::abs(func(q))*h/3 <= std::abs(result)*std::numeric_limits<double>::epsilon()){
            break;
        }
    }
    return result;
}

