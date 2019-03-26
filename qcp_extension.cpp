#include "qcp_extension.h"
#include "qcustomplot.h"
#include <cmath>
#include <limits>


static std::vector<QPen> color{QPen(Qt::blue), QPen(Qt::red), QPen(Qt::green),
                               QPen(Qt::darkYellow), QPen(Qt::darkMagenta), QPen(Qt::cyan)};

QCPRange get_x_range(QCustomPlot* plot){
    double mini = +std::numeric_limits<double>::infinity();
    double maxi = -std::numeric_limits<double>::infinity();

    for(int i=0; i<plot->plottableCount(); i++){
        bool found_range;
        auto p = plot->plottable(i);
        QCPRange range = p->getKeyRange(found_range);

        if(found_range){
            mini = std::min(mini, range.lower);
            maxi = std::max(maxi, range.upper);
        }
    }

    return QCPRange(mini, maxi);
}


void rescale_axes(QCustomPlot* plot, double x_lower, double x_upper, bool logplot){
    double mini = +std::numeric_limits<double>::infinity();
    double maxi = -std::numeric_limits<double>::infinity();

    bool found_any = false;

    auto key_range = QCPRange(x_lower, x_upper);
    for(int i=0; i<plot->plottableCount(); i++){
        bool found_range;
        auto p = plot->plottable(i);
        QCPRange range;
        if(logplot){
            range = p->getValueRange(found_range, QCP::sdPositive, key_range);
        } else {
            range = p->getValueRange(found_range, QCP::sdBoth, key_range);
        }
        if(found_range){
            mini = std::min(mini, range.lower);
            maxi = std::max(maxi, range.upper);
        }
        found_any = found_range || found_any;
    }

    plot->xAxis->setRange(key_range);
    if(found_any){
        double margin = 0.2;
        double relative_margin;
        if(logplot){
            double diff = maxi/mini;
            relative_margin = std::pow(diff,margin);
            plot->yAxis->setRange(mini/relative_margin, maxi*relative_margin);

            QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
            plot->yAxis->setTicker(logTicker);
        } else {
            double diff = maxi-mini;
            relative_margin = diff*margin;
            plot->yAxis->setRange(mini-relative_margin, maxi+relative_margin);
        }
    }

    plot->replot();
}


void add_graph(QCustomPlot* plot, std::vector<double> xvals, std::vector<double> vals,
              bool logplot, int parts){
    QVector<double> x = QVector<double>::fromStdVector(xvals);
    QVector<double> y = QVector<double>::fromStdVector(vals);

    if(logplot){
        plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    } else {
        plot->yAxis->setScaleType(QCPAxis::stLinear);
    }

    int graph_num = plot->plottableCount();
    int limit = 6*parts;
    if(graph_num >= 2*limit){
        for(int i=1; i<=limit; i++){
            plot->removePlottable(graph_num-i);
        }
        graph_num -= limit;
    }

    if (graph_num >= limit){
        graph_num -= limit;
    }

    QCPCurve* curve = new QCPCurve(plot->xAxis, plot->yAxis);

    auto pen = color[size_t(graph_num/parts)];
    pen.setStyle(Qt::SolidLine);
    curve->setPen(pen);
    if(parts>1 && graph_num%2 == 0){
        auto pen = curve->pen();
        pen.setStyle(Qt::DotLine);
        curve->setPen(pen);
    }
    curve->setData(x,y);

    if(logplot){
        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        plot->yAxis->setTicker(logTicker);
    }

    plot->replot();
}
