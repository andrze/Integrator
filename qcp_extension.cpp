#include "qcp_extension.h"
#include "qcustomplot.h"
#include <cmath>
#include <limits>


static std::vector<QPen> color{QPen(Qt::blue), QPen(Qt::red), QPen(Qt::darkGreen),
                               QPen(Qt::darkYellow), QPen(Qt::darkMagenta), QPen(Qt::darkCyan)};

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

std::pair<double, double> get_val_range(QCustomPlot* plot, double x_lower, double x_upper, bool logplot){

    double mini = +std::numeric_limits<double>::infinity();
    double maxi = -std::numeric_limits<double>::infinity();

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
    }

    return std::make_pair(mini,maxi);
}


void rescale_axes(QCustomPlot* plot, double x_lower, double x_upper, bool logplot){
    auto range = get_val_range(plot, x_lower, x_upper, logplot);

    double mini = range.first;
    double maxi = range.second;

    bool found_any = (mini < std::numeric_limits<double>::infinity()
            && maxi >-std::numeric_limits<double>::infinity());
    plot->xAxis->setRange(QCPRange(x_lower,x_upper));
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

void rescale_axes(QCustomPlot* plot, std::pair<double,double> x_range, std::pair<double,double> y_range, bool logplot){
    plot->xAxis->setRange(x_range.first,x_range.second);
    plot->yAxis->setRange(y_range.first,y_range.second);

    if(logplot){
        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        plot->xAxis->setTicker(logTicker);
        plot->yAxis->setTicker(logTicker);
    }
    plot->replot();
}

void add_graph(QCustomPlot* plot, std::vector<double> xvals, std::vector<double> vals,
              bool logplot, int parts){
    QVector<double> x = QVector<double>(xvals.begin(),xvals.end());
    QVector<double> y = QVector<double>(vals.begin(),vals.end());

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
    if(parts>1 && graph_num%2 == 1){
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


std::pair<std::vector<double>, std::vector<double> > posneg_part(std::vector<double> vals){

    std::vector<double> positive, negative;

    for(auto && v: vals){

        if(v > 0){
            positive.push_back(v);
            negative.push_back(0.);
        } else {
            positive.push_back(0);
            negative.push_back(-v);
        }
    }
    return std::make_pair(positive, negative);
}


