#include "realvector.h"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <physics.h>

RealVector::RealVector()
{
}

RealVector::RealVector(std::vector<double> coords){
    this->coords = coords;
}

RealVector& RealVector::operator+= (RealVector rhs) {
    if(coords.size() != rhs.coords.size()){
        throw std::invalid_argument( "Received RealVectors of different dimensions" );
    }

    for(size_t i=0; i<coords.size(); i++){
        coords[i] += rhs[i];
    }

    return *this;
}

RealVector& RealVector::operator*= (double rhs) {
    for(size_t i=0; i<coords.size(); i++){
        coords[i] *= rhs;
    }
    return *this;
}

RealVector& RealVector::operator-= (RealVector rhs) {
    return (*this) += ( (-1) * rhs );
}

RealVector& RealVector::operator/= (double rhs) {
    return (*this) *= (1/rhs);
}

double& RealVector::operator[] (size_t i){
    return coords[i];
}

RealVector operator + (RealVector lhs, RealVector rhs) {
    return lhs += rhs;
}

RealVector operator - (RealVector lhs, RealVector rhs) {
    return lhs -= rhs;
}

RealVector operator * (double lhs, RealVector rhs) {
    return rhs *= lhs;
}

RealVector operator * (RealVector lhs, double rhs) {
    return lhs *= rhs;
}

RealVector operator / (RealVector lhs, double rhs) {
    return lhs /= rhs;
}

double operator * (RealVector lhs, RealVector rhs) {
    if(lhs.coords.size() != rhs.coords.size()){
        throw std::invalid_argument( "Received RealVectors of different dimensions" );
    }

    double dot_prod = 0.;
    for(size_t i=0; i<lhs.coords.size(); i++){
        dot_prod += lhs[i] * rhs[i];
    }

    return dot_prod;
}

std::ostream& operator << (std::ostream& out, RealVector v) {
    for(size_t i=0; i<v.coords.size(); i++){
        out << v[i] << " ";
    }
    return out;
}

EquationSet::EquationSet(){

}

EquationSet::EquationSet(std::vector<std::function<double (RealVector, double, double)> > equations,
                         std::vector<std::function<double (double, double)> >  scale){
    this->equations = equations;
    this->scale = scale;
}

RealVector EquationSet::evaluate(RealVector point, double d){

    if(point.coords.size() != equations.size()){
        throw std::invalid_argument( "Equations and point of evaluation have different dimensions" );
    }

    std::vector<double> coords;
    double e = eta(point);

    for(size_t i=0; i<point.coords.size(); i++){
        coords.push_back( equations[i](point, e, d) );
    }

    return RealVector(coords);
}

RealVector filter(RealVector v, std::vector<bool> filter_vec){
    if(v.coords.size() != filter_vec.size()){
        throw std::invalid_argument( "Point and filter have different dimensions" );
    }

    for(size_t i=0; i<v.coords.size(); i++){
        if(!filter_vec[i]){
            v[i] = 0.;
        }
    }
    return v;
}

