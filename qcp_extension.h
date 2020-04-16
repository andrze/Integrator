#ifndef QCP_EXTENSION_H
#define QCP_EXTENSION_H
#include "qcustomplot.h"

QCPRange get_x_range(QCustomPlot* plot);

std::pair<double, double> get_val_range(QCustomPlot* plot, double x_lower, double x_upper, bool logplot=false);

void rescale_axes(QCustomPlot* plot, double x_lower, double x_upper, bool logplot=false);

void rescale_axes(QCustomPlot* plot, std::pair<double,double> x_range, std::pair<double,double> y_range, bool logplot=false);

void add_graph(QCustomPlot* plot, std::vector<double> xvals, std::vector<double> vals,
              bool logplot=false, int parts=1);

std::pair<std::vector<double>, std::vector<double> > posneg_part(std::vector<double> vals);

#endif // QCP_EXTENSION_H
