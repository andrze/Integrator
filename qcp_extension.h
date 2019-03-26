#ifndef QCP_EXTENSION_H
#define QCP_EXTENSION_H
#include "qcustomplot.h"

QCPRange get_x_range(QCustomPlot* plot);

void rescale_axes(QCustomPlot* plot, double x_lower, double x_upper, bool logplot=false);

void add_graph(QCustomPlot* plot, std::vector<double> xvals, std::vector<double> vals,
              bool logplot=false, int parts=1);


#endif // QCP_EXTENSION_H
