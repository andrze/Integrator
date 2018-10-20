#include<cmath>
#include<iostream>
#include<vector>
#include<functional>

double integral(std::function<double(double)> func){
    double h0=0.001;
    double q=0.;
    double result=0.;
    int num_steps=51;

    for(int i=0; i<100; i++){
        double h = h0/(std::pow(func(q),0.2));
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

        if(result+func(q)*h/3 == result){
            break;
        }
    }
    return result;
}

