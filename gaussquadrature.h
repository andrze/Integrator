#ifndef GAUSSQUADRATURE_H
#define GAUSSQUADRATURE_H

#include <vector>
#include <array>
#include <functional>

class GaussQuadrature
{
public:
    GaussQuadrature();
    double integrate(double a, double b, std::function<double(double)> f, size_t n);

private:
    std::array<std::vector<double>, 21> roots;
    std::array<std::vector<double>, 21> weights;
};

double gauss_legendre_integrate(std::function<double(double)> f);

#endif // GAUSSQUADRATURE_H
