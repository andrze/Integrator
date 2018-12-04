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

    int graph_num = plot->graphCount();
    int limit = 6*parts;
    if(graph_num >= 2*limit){
        for(int i=1; i<=limit; i++){
            plot->removeGraph(graph_num-i);
        }
        graph_num -= limit;
    }

    if (graph_num >= limit){
        graph_num -= limit;
    }

    plot->addGraph();
    plot->graph(graph_num)->setPen(color[size_t(graph_num/parts)]);
    if(parts>1 && graph_num%2 == 0){
        auto pen = plot->graph(graph_num)->pen();
        pen.setStyle(Qt::DotLine);
        plot->graph(graph_num)->setPen(pen);
    }
    plot->graph(graph_num)->setData(x,y);

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
            plot->yAxis->setRange(0,5);
        } else {
            double diff = max / min;
            if(diff <= 1.){
                diff = exp(5.);
            }
            diff = std::pow(diff, 0.1);
            plot->yAxis->setRange(std::min(min/diff, yrange.lower), std::max(max*diff, yrange.upper));
        }
    }
    plot->replot();
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    plots = std::vector<QCustomPlot*>{ui->alphaPlot, ui->msPlot, ui->mpPlot, ui->zPlot, ui->universalPlot,
                ui->uPlot, ui->tauPlot, ui->yPlot, ui->stiffPlot, ui->etaPlot,};

    EquationSet equations(ui->cubicRadioButton->isChecked());
    integrator = Integrator(equations);
    ui->tsenseBox->valueChanged(ui->tsenseBox->value());
    this->reset_xAxis();
    connect(this, SIGNAL(calculated(std::vector<Plot>*)), this, SLOT(plot_result(std::vector<Plot>*)));
}

MainWindow::~MainWindow()
{
    delete ui;
}

std::vector<Plot> MainWindow::integrate(){
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

std::vector<Plot> MainWindow::integrate(RealVector start_point){

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
            if(k < results[i][0].times.size()){
                file << results[i][0].times[k] << ',';
                for(size_t j=0; j<results[0].size(); j++){

                    file << results[i][j].vals[k] << ',';
                }
                file << (-exp(-results[i][0].times[k])*results[i][4].derivatives[k]/results[i][4].vals[k])<< ',';
            } else {
                for(size_t j=0; j<results[0].size()+2; j++){
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
        std::vector<Plot>* last_el;
        std::vector<Plot> p = integrate();
        {
            std::lock_guard<std::mutex> guard(this->data_mutex);
            results.push_back(p);
            last_el = &results.back();
        }

        std::string phase;

        int diag = phase_diagnosis(p);
        if(diag == 2){
            phase = "Uporządkowana";
        } else if(diag == 1){
            phase = "KT";
        } else {
            phase = "Nieuporządkowana";
        }


        std::cout<<"Obliczono, końcowa faza: " << phase <<std::endl;
        emit calculated(last_el);
    };

    std::thread t(func);

    t.detach();

}

void MainWindow::plot_result(std::vector<Plot>* res){
    auto p = *res;
    std::array<std::vector<double>, 7> plot_vals;

    addGraph(ui->alphaPlot, p[0].times, p[0].vals, true);

    addGraph(ui->msPlot, p[1].times, p[1].vals, true);

    addGraph(ui->mpPlot, p[2].times, p[2].vals, true);

    addGraph(ui->zPlot, p[3].times, p[3].vals, true, 2);
    addGraph(ui->zPlot, p[4].times, p[4].vals, true, 2);


    for(size_t i=0; i<p[0].vals.size(); i++){
        double a2 = std::pow(p[0].vals[i],2);
        plot_vals[0].push_back(p[1].vals[i]/a2);
        plot_vals[1].push_back(p[2].vals[i]/a2);
        plot_vals[2].push_back((p[3].vals[i]-p[4].vals[i])/a2);
        plot_vals[3].push_back(a2*p[4].vals[i]);
        plot_vals[4].push_back(-exp(-p[4].times[i])*p[4].derivatives[i]/p[4].vals[i]);
        plot_vals[5].push_back(-exp(-p[4].times[i])*p[4].derivatives[i]*a2/p[5].vals[i]);
        plot_vals[6].push_back(p[1].vals[i]/a2*exp(2*p[4].times[i])/std::pow(p[3].vals[i],2));
    }

    std::vector<QCustomPlot*> composite_plots{ui->uPlot, ui->tauPlot, ui->yPlot, ui->stiffPlot, ui->etaPlot, ui->universalPlot};
    for(size_t i=0; i<composite_plots.size(); i++){
        composite_plots[i]->yAxis->setScaleType(QCPAxis::stLinear);
        addGraph(composite_plots[i], p[0].times, plot_vals[i]);
    }

    //auto max = std::max_element(plot_vals[3].begin(), plot_vals[3].end());

    ui->etaAlphaPlot->xAxis->setRange(0,1.);
    ui->uAlphaPlot->xAxis->setRange(0,1.);
    addGraph(ui->etaAlphaPlot, plot_vals[3], plot_vals[4]);
    addGraph(ui->uAlphaPlot, plot_vals[3], plot_vals[6], true);

}

void MainWindow::on_clearButton_clicked()
{
    std::vector<QCustomPlot*> plots2 = plots;
    plots2.push_back(ui->etaAlphaPlot);
    plots2.push_back(ui->uAlphaPlot);
    for(auto&& p: plots2){
        p->clearGraphs();
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

void MainWindow::on_cubicRadioButton_clicked(){
    integrator = Integrator(EquationSet(ui->cubicRadioButton->isChecked()));
}

void MainWindow::on_hexRadioButton_clicked(){
    integrator = Integrator(EquationSet(ui->cubicRadioButton->isChecked()));
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
        int phase = phase_diagnosis(p);
        if(ui->orderedRadioButton->isChecked()){
            go_to_lower = (phase <= 1);
        } else {
            go_to_lower = (phase == 0);
        }

        if(go_to_lower){
            upper = T;
            T -= diff;
        } else {
            lower = T;
            T += diff;
        }
        ui->TSpinBox->setValue(T);

        std::cout<<ui->TSpinBox->value()<<std::endl;
    }
    double precision = 1e-6;

    ui->TSpinBox->setValue((upper+lower)/2);

    int i=0;
    while(std::abs(lower-upper) > precision){
        double T = ui->TSpinBox->value();

        std::cout<<i<<": "<<ui->TSpinBox->value()<<std::endl;

        auto p = this->integrate();

        bool go_to_lower;
        int phase = phase_diagnosis(p);
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

}


double MainWindow::find_fp(RealVector start_point){
    double diff = 0.002;

    std::cout<<std::setprecision(8);

    double lower = -1;
    double upper = -1;
    while(lower < 0 || upper < 0){
        auto p = this->integrate(start_point);

        bool go_to_lower;
        int phase = phase_diagnosis(p);
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

        std::cout<<start_point[5]<<std::endl;
    }
    double precision = 1e-7;

    start_point[5] = (upper+lower)/2;

    int i=0;
    while(std::abs(lower-upper) > precision){

        std::cout<<i<<": "<<start_point[5]<<' '<<lower<<' '<<upper<<std::endl;

        auto p = this->integrate(start_point);

        bool go_to_lower;
        int phase = phase_diagnosis(p);
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

    if(end <= 0. || start <= 0.){
        return;
    }
    size_t threads = 4;

    std::vector<double> points;
    for(size_t i=0; i<steps; i++){
        points.push_back(std::pow(end/start,double(i)/(steps-1)) * start);
    }

    auto func = [=](size_t first, size_t last, std::atomic<bool>* stop_ptr, std::promise<std::vector<std::pair<double,double> > > promises){
        std::vector<std::pair<double,double> > temperatures;
        for(size_t i=first; i<last; i++){
            double lambda = points[i];
            RealVector point = start_point;
            point[2] = lambda;
            double T = this->find_fp(point);
            std::cout<<"Punkt "<<i<<", Lambda: "<< lambda<<", Temparatura :"<<T<<std::endl;
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
