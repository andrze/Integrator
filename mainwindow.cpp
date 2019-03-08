#include <vector>
#include <functional>
#include <fstream>
#include <QString>
#include <QFileDialog>
#include <algorithm>
#include <cmath>
#include <array>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "integrator.h"
#include "realvector.h"
#include "physics_cubic.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <thread>

static std::vector<QPen> color{QPen(Qt::blue), QPen(Qt::red), QPen(Qt::green),
                               QPen(Qt::darkYellow), QPen(Qt::darkMagenta), QPen(Qt::cyan)};


void addGraph(QCustomPlot* plot, std::vector<double> xvals, std::vector<double> vals,
              bool logplot=false, int parts=1){
    QVector<double> x = QVector<double>::fromStdVector(xvals);
    QVector<double> y = QVector<double>::fromStdVector(vals);

    if(plot->yAxis->range().lower == 0. || plot->yAxis->range().lower == 5.){
        if(logplot){
            plot->yAxis->setRange(0.9, 1.1);
        } else {
            plot->yAxis->setRange(0., 0.00000001);
        }
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

    auto yrange = plot->yAxis->range();
    auto max = *std::max_element(vals.begin(), vals.end());
    if(!logplot){
        auto min = *std::min_element(vals.begin(), vals.end());
        double diff = max - min;
        if(diff > 0){
            plot->yAxis->setRange(std::min(min-diff/8, yrange.lower), std::max(max+diff/8, yrange.upper));
        }

    } else {
        auto min = max;
        for(auto v: vals){
            if(v < min && v > 0.){
                min = v;
            }
        }
        plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
        QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
        plot->yAxis->setTicker(logTicker);

        if(max <= 0){
            plot->yAxis->setRange(1e-10,10);
        } else {
            double diff = max / min;
            if(diff <= 1.){
                diff = exp(5.);
            }
            diff = std::pow(diff, 0.1);
            plot->yAxis->setRange(std::min(std::abs(min/diff), std::abs(yrange.lower)), std::max(max*diff, yrange.upper));
        }
    }
    plot->replot();
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    plots = std::vector<QCustomPlot*>{ui->alphaPlot, ui->mPlot, ui->kappaPlot, ui->zPlot, ui->universalPlot,
                ui->uPlot, ui->yPlot, ui->m2Plot, ui->etaPlot,};

    EquationSet equations(ui->cubicRadioButton->isChecked(), ui->dimensionBox->value());
    integrator = Integrator(equations);
    ui->tsenseBox->valueChanged(ui->tsenseBox->value());
    this->reset_xAxis();
    connect(this, SIGNAL(calculated(PlotSet*)), this, SLOT(plot_result(PlotSet*)));
}

MainWindow::~MainWindow()
{
    delete ui;
}

std::pair<PlotSet, int> MainWindow::integrate(){
    std::vector<double> coords({1,
                                ui->uBox->value(), ui->tauBox->value(),
                                1, 1,
                                ui->TSpinBox->value()});
    RealVector start(coords);

    return integrator.integrate(exp(-ui->startTimeBox->value()),
                                exp(-ui->endTimeBox->value()),
                                ui->deltaTBox->value(),
                                start);
}

std::pair<PlotSet, int> MainWindow::integrate(RealVector start_point){

    return integrator.integrate(exp(-ui->startTimeBox->value()),
                                exp(-ui->endTimeBox->value()),
                                ui->deltaTBox->value(),
                                start_point);
}

void MainWindow::on_fileButton_clicked(){
    {
        std::string filename = QFileDialog::getSaveFileName(this, tr("Save file"), "", tr("CSV files (*.csv)")).toStdString();
        if( filename.size() == 0 ){
            return;
        }
        if( filename.size() < 4 || filename.substr(filename.size()-4,4) != ".csv" ){
            filename += ".csv";
        }
        ui->saveButton->setEnabled(true);
        ui->filenameEdit->setText(tr(&filename[0]));
    }
}

void MainWindow::on_saveButton_clicked(){
    if(results.empty()){
        return;
    }

    std::ofstream file(ui->filenameEdit->text().toStdString(), std::ios::out);
    file<<std::setprecision(10);
    //file<<std::setprecision(std::numeric_limits<long double>::digits10 + 1);

    std::vector<std::string> headers{"s", "Alpha", "Sigma Mass", "Pi Mass", "Z Sigma", "Z Pi", "T", "Eta"};

    for(size_t j=0; j<results.size(); j++){
        for(size_t i=0; i<headers.size(); i++){
            file << headers[i] << ',';
        }
        file << ',';
    }
    file << '\n';

    size_t max_length = results[0][0].times.size();
    for(size_t i=0; i<results.size(); i++){
        max_length = std::max(results[i][0].times.size(), max_length);
    }

    for(size_t k=0; k<max_length; k++){
        for(size_t i=0; i<results.size(); i++){
            if(k < results[i].plot_size()){
                file << results[i][0].times[k] << ',';
                for(size_t j=0; j<results[i].plot_number(); j++){

                    file << results[i][j].values[k] << ',';
                }
                file << (-exp(-results[i][4].times[k])*results[i][4].derivatives[k]/results[i][4].values[k])<< ',';
            } else {
                for(size_t j=0; j<results[i].plot_number()+2; j++){
                    file << ',';
                }
            }
            file << ',';
        }
       file << '\n';
    }
    file.close();

}

void MainWindow::on_tsenseBox_valueChanged(int arg1){
    this->ui->tauBox->setSingleStep(pow(10.,-arg1));
}

void MainWindow::on_pushButton_clicked(){

    auto func = [=](){
        PlotSet* last_el;
        ui->dimensionBox->setEnabled(false);
        auto p = integrate();
        {
            std::lock_guard<std::mutex> guard(this->data_mutex);
            results.push_back(p.first);
            last_el = &results.back();
        }

        std::string phase;

        int diag = p.second;
        if(diag == 2){
            phase = "Uporządkowana";
        } else if(diag == 1){
            phase = "KT";
        } else {
            phase = "Nieuporządkowana";
        }


        std::cout<<"Obliczono, końcowa faza: " << phase <<std::endl;
        emit calculated(last_el);

        ui->dimensionBox->setEnabled(true);
    };

    std::thread t(func);

    t.detach();

}

void MainWindow::plot_result(PlotSet *res){
    auto p = *res;
    std::array<std::vector<double>, 9> plot_vals;
    std::vector<double> time, log_time;
    for(auto t: p[0].times){
        time.push_back(t);
        log_time.push_back(-log(t));
    }
    double d = ui->dimensionBox->value();

    if(p[2].values[0] != 0){
        addGraph(ui->mPlot, log_time, p[1].values, true, 2);
        addGraph(ui->mPlot, log_time, p[2].values, true, 2);
    } else {
        addGraph(ui->mPlot, log_time, p[1].values, true, 1);
    }
    addGraph(ui->zPlot, log_time, p[3].values, true, 2);
    addGraph(ui->zPlot, log_time, p[4].values, true, 2);


    for(size_t i=0; i<p[0].values.size(); i++){
        double pow_time = 1;
        if(d!=2){
            pow_time = std::pow(time[i],(2-d)/2);
        }
        double a2 = std::pow(p[0].values[i],2);
        double kappa = a2*p[4].values[i]*pow_time*pow_time;
        plot_vals[0].push_back((p[3].values[i]-p[4].values[i])/kappa);
        plot_vals[1].push_back(kappa);
        plot_vals[2].push_back(p.eta(i));
        plot_vals[3].push_back(p.eta(i)*kappa/p[5].values[i]);
        plot_vals[4].push_back(a2);
        plot_vals[5].push_back(p[1].values[i]/a2 * std::pow(time[i], d-4) * std::pow(p[4].values[i],-2));
        plot_vals[6].push_back(p[2].values[i]/a2 * std::pow(time[i], d-4) * std::pow(p[4].values[i],-2));
        plot_vals[7].push_back(p[1].values[i]/p[4].values[i]/time[i]/time[i]);
        plot_vals[8].push_back(p[2].values[i]/p[4].values[i]/time[i]/time[i]);
    }

    std::vector<QCustomPlot*> composite_plots{ui->yPlot, ui->kappaPlot, ui->etaPlot, ui->universalPlot, ui->alphaPlot};
    std::vector<bool> logplot{                true,      true,         false,       false,             true};
    for(size_t i=0; i<composite_plots.size(); i++){
        if(logplot[i]){
            composite_plots[i]->yAxis->setScaleType(QCPAxis::stLogarithmic);
        } else {
            composite_plots[i]->yAxis->setScaleType(QCPAxis::stLinear);
        }
        addGraph(composite_plots[i], log_time, plot_vals[i], logplot[i]);
    }

    ui->m2Plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    addGraph(ui->m2Plot, log_time, plot_vals[7], true, 2);
    addGraph(ui->m2Plot, log_time, plot_vals[8], true, 2);

    ui->uPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    if(plot_vals[6][0] != 0){
        addGraph(ui->uPlot, log_time, plot_vals[5], true, 2);
        addGraph(ui->uPlot, log_time, plot_vals[6], true, 2);
    } else {
        addGraph(ui->uPlot, log_time, plot_vals[5], true, 1);
    }
    //auto max = std::max_element(plot_vals[3].begin(), plot_vals[3].end());

    ui->etaAlphaPlot->xAxis->setRange(0,1.);
    ui->uAlphaPlot->xAxis->setRange(0,1.);
    addGraph(ui->etaAlphaPlot, plot_vals[1], plot_vals[2]);
    addGraph(ui->uAlphaPlot, plot_vals[1], plot_vals[0], true);


    double lower=ui->ulPlot->xAxis->range().lower;
    double upper=ui->ulPlot->xAxis->range().upper;
    if(lower == 0){
        lower = INFINITY;
        upper = -INFINITY;
    }
    lower = std::min(lower, *std::min_element(plot_vals[5].begin(), plot_vals[5].end()));
    upper = std::max(upper, *std::max_element(plot_vals[5].begin(), plot_vals[5].end()));

    ui->ulPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    ui->ulPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);

    ui->ulPlot->xAxis->setRange(lower, upper);
    QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
    ui->ulPlot->yAxis->setTicker(logTicker);
    ui->ulPlot->xAxis->setTicker(logTicker);
    addGraph(ui->ulPlot, plot_vals[5], plot_vals[6], true);

    std::cout<<"\a"<<std::flush;
}

void MainWindow::on_clearButton_clicked()
{
    std::vector<QCustomPlot*> plots2 = plots;
    plots2.push_back(ui->etaAlphaPlot);
    plots2.push_back(ui->uAlphaPlot);
    plots2.push_back(ui->ulPlot);
    for(auto&& p: plots2){
        p->clearPlottables();
        p->replot();
    }
    results.clear();
}

void MainWindow::reset_xAxis(){

    for(auto&& p: plots){
        double lower = ui->xMinSpinBox->value();
        double upper = ui->xMaxSpinBox->value();
        p->xAxis->setRange(lower, upper);
        p->replot();
    }
}

void MainWindow::on_xMinSpinBox_valueChanged(double){
    reset_xAxis();
}

void MainWindow::on_xMaxSpinBox_valueChanged(double){
    reset_xAxis();
}

void MainWindow::restart_integrator(){
    EquationSet equations(ui->cubicRadioButton->isChecked(), ui->dimensionBox->value());
    integrator = Integrator(equations);
}

void MainWindow::on_cubicRadioButton_clicked(){
    restart_integrator();
}

void MainWindow::on_hexRadioButton_clicked(){
    restart_integrator();
}

void MainWindow::on_fpButton_clicked(){
    double diff = 0.002;

    std::cout<<std::setprecision(8);

    double lower = -1;
    double upper = -1;
    while(lower < 0 || upper < 0){
        double T = ui->TSpinBox->value();
        auto p = this->integrate();

        bool go_to_lower;
        int phase = p.second;
        if(ui->orderedRadioButton->isChecked()){
            go_to_lower = (phase <= 1);
        } else {
            go_to_lower = (phase == 0);
        }

        if(go_to_lower){
            upper = T;
            T -= diff;
            if(T > 4.1*diff){
                diff *=2;
            }
            while(T <= diff){
                diff /= 2;
            }
        } else {
            lower = T;
            T += diff;
            diff *= 2;
        }
        ui->TSpinBox->setValue(T);

        std::cout<<ui->TSpinBox->value()<<std::endl;
    }
    double precision = std::pow(10,-ui->tsenseBox->value());

    ui->TSpinBox->setValue((upper+lower)/2);

    int i=0;
    while(std::abs(lower-upper) > precision){
        double T = ui->TSpinBox->value();

        std::cout<<i<<": "<<ui->TSpinBox->value()<<std::endl;

        auto p = this->integrate();

        bool go_to_lower;
        int phase = p.second;
        if(ui->orderedRadioButton->isChecked()){
            go_to_lower = (phase <= 1);
        } else {
            go_to_lower = (phase == 0);
        }


        if(go_to_lower){
            upper = T;
        } else {
            lower = T;
        }
        ui->TSpinBox->setValue((lower + upper)/2);
        i++;
    }

    std::cout<<"\a"<<std::flush;
}


double MainWindow::find_fp(RealVector start_point){
    double diff = 0.002;

    std::cout<<std::setprecision(8);

    double lower = -1;
    double upper = -1;
    size_t steps = 0;
    while(lower < 0 || upper < 0){
        auto p = this->integrate(start_point);

        int phase = p.second;

        bool go_to_lower;
        if(ui->orderedRadioButton->isChecked()){
            go_to_lower = (phase <= 1);
        } else {
            go_to_lower = (phase == 0);
        }

        if(go_to_lower){
            upper = start_point[5];
            start_point[5] -= diff;
        } else {
            lower = start_point[5];
            start_point[5] += diff;
        }

        std::cout<<"T: "<<start_point[5]<<" L: "<<start_point[2]<<std::endl;
        steps++;
        if(steps>50){
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    double precision = 1e-8;

    start_point[5] = (upper+lower)/2;

    int i=0;
    while(std::abs(lower-upper) > precision){

        std::cout<<i<<" T: "<<start_point[5]<<" L: "<<start_point[2]<<std::endl;

        auto p = this->integrate(start_point);

        bool go_to_lower;
        int phase = p.second;
        if(ui->orderedRadioButton->isChecked()){
            go_to_lower = (phase <= 1);
        } else {
            go_to_lower = (phase == 0);
        }


        if(go_to_lower){
            upper = start_point[5];
        } else {
            lower = start_point[5];
        }
        start_point[5] = (lower + upper)/2;
        i++;
    }
    return start_point[5];
}

void MainWindow::on_CLButton_clicked()
{
    double start = ui->tauBox->value();
    double end = ui->criticalLineSpinBox->value();
    size_t steps = size_t(ui->criticalNumberBox->value())+1;
    RealVector start_point{std::vector<double>{1,ui->uBox->value(),0,1,1,ui->TSpinBox->value()}};

    bool logarithmic = ui->hexRadioButton->isChecked();

    if(end <= 0. || start <= 0.){
        return;
    }
    size_t threads = 4;

    std::vector<double> points;

    for(size_t i=0; i<steps; i++){
        if(logarithmic){
            points.push_back(std::pow(end/start,double(i)/(steps-1)) * start);
        } else {
            points.push_back((end-start)*double(i)/(steps-1)+start);
        }
    }

    auto func = [=](size_t first, size_t last, std::atomic<bool>* stop_ptr, std::promise<std::vector<std::pair<double,double> > > promises){
        std::vector<std::pair<double,double> > temperatures;
        double T = start_point.coords[5];
        for(size_t i=first; i<last; i++){
            double lambda = points[i];
            RealVector point = start_point;
            point[2] = lambda;
            point[5] = T;
            double T = this->find_fp(point);
            std::cout<<"Punkt "<<i<<", Lambda: "<< lambda<<", Temperatura :"<<T<<std::endl;
            temperatures.push_back(std::make_pair(T, lambda));
            if(*stop_ptr){
                *stop_ptr = false;
                break;
            }
        }
        promises.set_value(temperatures);
    };

    std::vector<std::thread > calculations;
    std::vector<std::promise<std::vector<std::pair<double, double> > > > promises;
    std::vector<std::future<std::vector<std::pair<double, double> > > > values;

    for(size_t i=0; i<threads; i++){
        std::promise<std::vector<std::pair<double, double> > > promise;
        values.push_back(promise.get_future());
        calculations.push_back(std::thread(func, (steps*i)/threads, (steps*(i+1))/threads, &stop,
                                           std::move(promise)));
    }

    std::vector<std::pair<double, double> > temperatures;

    for(size_t i=0; i<threads; i++){
        values[i].wait();
        calculations[i].join();
        std::vector<std::pair<double, double> > result = values[i].get();
        temperatures.insert(temperatures.end(), result.begin(), result.end());
    }

    std::ofstream file("/Users/andrzejchlebicki/Desktop/critline.csv", std::ios::app);
    file << std::setprecision(8);
    file << "T, lambda\n";
    for(auto && l : temperatures){
        file << l.first << ',' << l.second <<std::endl;
    }

    return;
}

void MainWindow::on_stopButton_clicked(){
    this->stop = true;
}

void MainWindow::on_dimensionBox_valueChanged(double)
{
    restart_integrator();
}
