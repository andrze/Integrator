#ifndef PLOTSET_H_
#define PLOTSET_H_

#include "realvector.h"
#include "plot.h"
#include "equationset.h"

class PlotSet {
public:
    PlotSet();
    PlotSet(size_t size, double d);
    virtual ~PlotSet();

    void push_to_each(RealVector values, RealVector ders, double t);
    void pop_from_each();
    RealVector point(size_t i, bool values=true);
    RealVector starting_point();
    RealVector back_vals();
    RealVector back_ders();
    double back_time();

    size_t plot_size();
    Plot& operator[](size_t i);
    size_t plot_number();

    double eta(size_t k);
    std::vector<double> eta_plot();
    std::pair<double,double> rescaled(size_t plot, size_t pos);
    Plot rescaled(size_t k);
    int phase_diagnosis();
    int phase = -1;

private:
    std::vector<Plot> plots;
    double d;

};

std::ostream& operator<<(std::ostream& out, PlotSet plots);

#endif /* PLOTSET_H_ */
