/*
 * PlotSet.cpp
 *
 *  Created on: 15 gru 2018
 *      Author: andrzejchlebicki
 */

#include "plotset.h"
#include <limits>
#include <stdexcept>
#include <cmath>
#include <iostream>

PlotSet::PlotSet(size_t size, double d) {
    for(size_t j=0; j<size; j++){
		plots.push_back(Plot());
	}
    this->d = d;
}

PlotSet::~PlotSet() {
}

void PlotSet::push_to_each(RealVector values, RealVector ders, double t){
	for(size_t i=0; i<plots.size(); i++){
		plots[i].values.push_back(values[i]);
		plots[i].derivatives.push_back(ders[i]);
		plots[i].times.push_back(t);
	}
}

void PlotSet::pop_from_each(){
	for(size_t i=0; i<plots.size(); i++){
		plots[i].values.pop_back();
		plots[i].derivatives.pop_back();
		plots[i].times.pop_back();
	}
}

RealVector PlotSet::point(size_t i, bool values){
    std::vector<double> coords;
    for(size_t j=0; j<plots.size(); j++){
        if(values){
            coords.push_back(plots[j].values[i]);
        } else {
            coords.push_back(plots[j].derivatives[i]);
        }
    }
    return RealVector(coords);
}


RealVector PlotSet::starting_point(){
    return point(0);
}

RealVector PlotSet::back_vals(){
    return point(plot_size()-1, true);
}

RealVector PlotSet::back_ders(){
    return point(plot_size()-1, false);
}

double PlotSet::back_time(){
	double t=plots[0].times.back();
	double precision = t*std::numeric_limits<double>::epsilon();
	for(size_t i=1; i<plots.size(); i++){
		if(std::abs(t - plots[i].times.back()) > precision){
			throw std::invalid_argument("Simulation time different in plots on same level");
		}
	}
	return t;
}

size_t PlotSet::plot_size(){
	size_t size = plots[0].size();
	for(size_t i=1; i<plots.size(); i++){
		if(plots[i].size() != size){
			throw std::invalid_argument("Sizes of plots in PlotSet differ");
		}
	}
	return size;
}

Plot& PlotSet::operator[](size_t i){

	return plots[i];
}

size_t PlotSet::plot_number(){
	return plots.size();
}

double PlotSet::eta(size_t k){
    return plots[5].log_der(k);
}

std::pair<double,double> PlotSet::rescaled(size_t plot, size_t pos){
    if(plot > this->plot_number()){
        throw std::invalid_argument("Argument larger than number of plots");
    }

    if(this->plot_number() < 5){
        throw std::invalid_argument("Not enough plots to allow rescaling");
    }

    if(pos >= plots[plot].size()){
        throw std::invalid_argument("Rescaled position out of plot scope");
    }

    double Z = plots[5].values[pos];
    double eta = plots[5].log_der(pos);

    Plot rescaled_plot = plots[plot];
    double time = rescaled_plot.times[pos];
    double scaling=1;
    double scaling_der=0;

    if(plot==0){ // Rescaling for kappa
        scaling = Z*std::exp(time*(-d+2));
        scaling_der = scaling*(-d+2-eta);
    } else if(plot==1 || plot==2){ // Rescaling for u and lambda
        scaling = std::pow(Z,-2)*std::exp(time*(d-4));
        scaling_der = scaling*(d-4+2*eta);
    } else if(plot==3 || plot==6){ // Rescaling for v and J
        scaling = std::pow(Z,-3)*std::exp(time*(2*d-6));
        scaling_der = scaling*(2*d-6+3*eta);
    } else if(plot==4){  // Rescaling for Y
        scaling = std::pow(Z,-2)*std::exp(time*(d-2));
        scaling_der = scaling*(d-2+2*eta);
    } else { // Rescaling for Z
        scaling = 1/Z;
        scaling_der = -eta*scaling;
    }
    double val = rescaled_plot.values[pos]*scaling;
    double der = rescaled_plot.derivatives[pos]*scaling+rescaled_plot.values[pos]*scaling_der;
    return std::make_pair(val, der);
}

Plot PlotSet::rescaled(size_t plot){
    if(plot > this->plot_number()){
        throw std::invalid_argument("Argument larger than number of plots");
    }

    if(this->plot_number() < 5){
        throw std::invalid_argument("Not enough plots to allow rescaling");
    }


    Plot rescaled_plot;
    rescaled_plot.times = this->plots[plot].times;
    for(size_t i=0; i<this->plot_size(); i++){
        std::pair<double, double> rescaled_point = rescaled(plot, i);

        rescaled_plot.values.push_back(rescaled_point.first);
        rescaled_plot.derivatives.push_back(rescaled_point.second);
    }

    return rescaled_plot;
}

int PlotSet::phase_diagnosis(){
    size_t plot_size = this->plot_size();

    size_t start, end;
    if(plot_size < 21){
		start = 1;
    } else {
        start = plot_size-20;
    }
    if(plot_size < 11){
        end = 0;
    } else {
        end = plot_size-10;
    }

    for(size_t i=start; i<end; i++){
        bool ordered = true;

        for(size_t j=0; j<plot_number(); j++){
            auto rescaling = rescaled(j, i);
            if(j==5) continue;
            if(rescaling.first < std::numeric_limits<double>::epsilon()) continue;
            if(j>0 && std::abs(plots[2].values.front()) < std::numeric_limits<double>::epsilon()) continue;

            double log_der = std::abs(rescaling.second/rescaling.first);
            if(log_der>1e-04){
                ordered = false;
            }

        }
        if(ordered){
            return 2;
        }
    }
    return 0;
}

std::ostream& operator<<(std::ostream& out, PlotSet plots){
	for(size_t i=0; i<plots.plot_size(); i++){
		for(size_t j=0; j<plots.plot_number(); j++){
			auto slice = plots[j][i];
			if(j == 0){
				out << -log(slice[0]) << ", ";
			}
			out << slice[1] << ", " << slice[2] << ", ";
		}
		out << "\n";
	}
	return out;
}
