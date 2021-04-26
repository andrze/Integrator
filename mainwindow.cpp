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
#include "physics.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <thread>
#include "qcp_extension.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    plots = std::vector<QCustomPlot*>{ui->alphaPlot, ui->lambdaPlot, ui->kappaPlot, ui->g3Plot, ui->zPlot,
                ui->uPlot, ui->yPlot, ui->etaPlot,ui->rho4Plot,ui->u21Plot, ui->mPlot, ui->m2Plot};

    is_plot_log = std::vector<bool>{true, true, true, true, false, true, true, false, true, true, true, true};

    EquationSet equations(ui->dimensionBox->value(), ui->NspinBox->value());
    integrator = Integrator(equations);
    ui->tsenseBox->valueChanged(ui->tsenseBox->value());
    this->reset_axes();
    connect(this, SIGNAL(calculated(PlotSet*)), this, SLOT(plot_result(PlotSet*)));
}

MainWindow::~MainWindow()
{
    delete ui;
}

PlotSet MainWindow::integrate(){
    std::vector<double> coords({ui->TSpinBox->value(),ui->uBox->value(),1});
    RealVector start(coords);

    return integrator.integrate(ui->startTimeBox->value(),
                                ui->endTimeBox->value(),
                                ui->deltaTBox->value(),
                                start);
}

PlotSet MainWindow::integrate(RealVector start_point){

    return integrator.integrate(ui->startTimeBox->value(),
                                ui->endTimeBox->value(),
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

    std::vector<std::string> headers{"s", "Kappa", "u", "l", "v", "Y", "Z", "J", "Eta"};

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
            auto res = results[i];
            if(k < res.plot_size()){
                file << res[0].times[k] << ',';
                for(size_t j=0; j<res.plot_number(); j++){

                    file << res[j].values[k] << ',';
                }
                file << res[5].derivatives[i]/res[5].values[i]<< ',';
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
    this->ui->TSpinBox->setSingleStep(pow(10.,-arg1));
}

void MainWindow::on_pushButton_clicked(){
    auto func = [=](){
        PlotSet* last_el;
        ui->dimensionBox->setEnabled(false);
        auto p = integrate();
        {
            std::lock_guard<std::mutex> guard(this->data_mutex);
            results.push_back(p);
            last_el = &results.back();
        }

        std::string phase;

        int diag = p.phase;
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
    std::array<std::vector<double>, 11> plot_vals;
    std::vector<double> time = p[0].times;

    //add_graph(ui->vPlot, time, p[5].abs_values(), true, 2);
    //add_graph(ui->vPlot, time, p[6].abs_values(), true, 2);

    plot_vals[0] = p[0].values;           // kappa
    plot_vals[1] = p[2].values;           // Z
    plot_vals[2] = p.rescaled(0).values;  // alpha^2
    plot_vals[3] = p.eta_plot();
    plot_vals[7] = p.rescaled(1).abs_values();      // u
    //plot_vals[5] = p.rescaled(4).abs_values();      // lambda

    std::vector<double> u = p[1].abs_values();
    //std::vector<double> l = p[4].abs_values();
    for(size_t i=0; i<p[0].values.size(); i++){
        plot_vals[4].push_back(2*plot_vals[0][i]*u[i]);
        plot_vals[5].push_back(2*plot_vals[2][i]*plot_vals[7][i]);
        plot_vals[6].push_back(plot_vals[4][i]*u[i]);
    }

    //auto u40_split = posneg_part(p[7].values);
    //plot_vals[6] = u40_split.first;
    //plot_vals[7] = u40_split.second;

    //auto y_split = posneg_part(p[2].values);
    //plot_vals[8] = y_split.first;
    //plot_vals[9] = y_split.second;

    std::vector<QCustomPlot*> composite_plots{ui->kappaPlot, ui->zPlot, ui->alphaPlot, ui->etaPlot, ui->mPlot, ui->m2Plot, ui->g3Plot};
    std::vector<bool> logplot{                true,          true,             true, false, true, true, true};
    for(size_t i=0; i<composite_plots.size(); i++){
        add_graph(composite_plots[i], time, plot_vals[i], logplot[i]);
    }

    /*ui->m2Plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    add_graph(ui->m2Plot, p[5].abs_values(), p[6].abs_values(), true, 1);


    auto m2range = get_x_range(ui->m2Plot);
    m2range.upper = std::min(200., m2range.upper);
    rescale_axes(ui->m2Plot, m2range.lower, m2range.upper, true);
*/

    add_graph(ui->uPlot, time, p[1].abs_values(), true, 1);
    //add_graph(ui->uPlot, time, p[4].abs_values(), true, 2);

    add_graph(ui->lambdaPlot, time, plot_vals[7], true, 1);
    //add_graph(ui->lambdaPlot, time, plot_vals[5], true, 2);

    ui->yPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    add_graph(ui->yPlot, plot_vals[0], plot_vals[7], true);
    auto range = get_x_range(ui->yPlot);
    range.upper = std::min(200., range.upper);
/*
    add_graph(ui->rho4Plot, time, u40_split.first, true, 2);
    add_graph(ui->rho4Plot, time, u40_split.second, true, 2);

    add_graph(ui->u21Plot, time, p[8].abs_values(), true, 2);
    add_graph(ui->u21Plot, time, p[9].abs_values(), true, 2);

    ui->ulPlot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    add_graph(ui->ulPlot, p[3].abs_values(), p[4].abs_values(), true);
    auto range = get_x_range(ui->ulPlot);
    range.upper = std::min(200., range.upper);

    rescale_axes(ui->ulPlot, range.lower, range.upper, true);


    QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog),logTicker2(new QCPAxisTickerLog);
    ui->m2Plot->xAxis->setTicker(logTicker);
    ui->ulPlot->xAxis->setTicker(logTicker2);
*/
    reset_axes();
    std::cout<<"\a"<<std::flush;
}

void MainWindow::on_clearButton_clicked(){
    std::vector<QCustomPlot*> plots2 = plots;
    for(auto&& p: plots2){
        p->clearPlottables();
        p->replot();
    }
    results.clear();

}

void MainWindow::on_StabMatButton_clicked()
{
    double end_time = ui->endTimeBox->value();
    ui->endTimeBox->setValue(1);

    auto plots = integrate();

    RealVector base_der = plots.point(1, false);

    RealVector fp = plots.starting_point();

    std::vector<RealVector> stab_mat;
    std::vector<double> deviation = {0.001,0.001,-1.,0.001,0.001,0.001};
    for(size_t i=0; i<deviation.size(); i++){
        RealVector start = fp;
        double epsilon = start[i]*deviation[i];
        start[i] += epsilon;

        auto dev_plots = integrate(start);


        stab_mat.push_back((dev_plots.point(1, false)-base_der)/epsilon);
    }

    std::cout << "\nSTABILITY MATRIX\n\n";
    for(size_t i=0; i<deviation.size(); i++){
        std::cout<<stab_mat[i]<<'\n';
    }

    ui->endTimeBox->setValue(end_time);
}

void MainWindow::reset_axes(){

    for(size_t i=0; i<plots.size(); i++){
        double lower = ui->xMinSpinBox->value();
        double upper = ui->xMaxSpinBox->value();

        rescale_axes(plots[i], lower, upper, is_plot_log[i]);
    }
}

void MainWindow::on_xMinSpinBox_valueChanged(double){
    reset_axes();
}

void MainWindow::on_xMaxSpinBox_valueChanged(double){
    reset_axes();
}

void MainWindow::restart_integrator(){
    EquationSet equations(ui->dimensionBox->value(), ui->NspinBox->value());
    integrator = Integrator(equations);
}

void MainWindow::on_fpButton_clicked(){
    double diff = 0.002;

    std::cout<<std::setprecision(8);

    double lower = -1;
    double upper = -1;
    while(lower < 0 || upper < 0){
        double alpha = ui->TSpinBox->value();
        auto p = this->integrate();

        bool go_to_lower;
        int phase = p.phase;
        if(ui->orderedRadioButton->isChecked()){
            go_to_lower = (phase == 2);
        } else {
            go_to_lower = (phase != 1);
        }

        if(go_to_lower){
            upper = alpha;
            alpha -= diff;
            if(alpha > 4.1*diff){
                diff *=2;
            }
            while(alpha <= diff){
                diff /= 2;
            }
        } else {
            lower = alpha;
            alpha += diff;
            diff *= 2;
        }
        ui->TSpinBox->setValue(alpha);

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
        int phase = p.phase;
        if(ui->orderedRadioButton->isChecked()){
            go_to_lower = (phase == 2);
        } else {
            go_to_lower = (phase != 1);
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

        int phase = p.phase;

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
        int phase = p.phase;
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
/*
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
*/
void MainWindow::on_stopButton_clicked(){
    this->stop = true;
}

void MainWindow::on_dimensionBox_valueChanged(double)
{
    restart_integrator();
}


void MainWindow::on_NspinBox_valueChanged(double )
{
    restart_integrator();
}
