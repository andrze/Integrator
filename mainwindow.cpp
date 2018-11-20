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
#include <limits>

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
    if(graph_num >= limit){
        plot->clearGraphs();
        graph_num = 0;
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
            plot->yAxis->setRange(std::min(min-diff/4, yrange.lower), std::max(max+diff/4, yrange.upper));
        }

    } else {
        auto min = max;
        for(auto v: vals){
            if(v < min && v > 0.){
                min = v;
            }
        }
        plot->yAxis->setScaleType(QCPAxis::stLogarithmic);

        if(max <= 0){
            plot->yAxis->setRange(0,5);
        } else {
            double diff = max / min;
            if(diff <= 1.){
                diff = exp(5.);
            }
            diff = std::pow(diff, 0.2);
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

    EquationSet equations;
    integrator = Integrator(equations);
    ui->tsenseBox->valueChanged(ui->tsenseBox->value());
    this->reset_xAxis();
}

MainWindow::~MainWindow()
{
    delete ui;
}

std::vector<Plot> MainWindow::integrate(){
    std::vector<double> coords({1,
                                1, ui->tauBox->value(),
                                1, 1,
                                ui->TSpinBox->value()});
    RealVector start(coords);

    return integrator.integrate(exp(-ui->startTimeBox->value()),
                                exp(-ui->endTimeBox->value()),
                                ui->deltaTBox->value(),
                                start);
}

void MainWindow::on_fileButton_clicked()
{
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

void MainWindow::on_saveButton_clicked()
{

    std::ofstream file(ui->filenameEdit->text().toStdString());
    //file << results;
}

void MainWindow::on_tsenseBox_valueChanged(int arg1)
{
    this->ui->tauBox->setSingleStep(pow(10.,-arg1));
}

void MainWindow::on_pushButton_clicked(){

    std::vector<Plot> p = this->integrate();
    std::array<std::vector<double>, 6> plot_vals;

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
    }

    std::vector<QCustomPlot*> composite_plots{ui->uPlot, ui->tauPlot, ui->yPlot, ui->stiffPlot, ui->etaPlot, ui->universalPlot};
    for(size_t i=0; i<composite_plots.size(); i++){
        composite_plots[i]->yAxis->setScaleType(QCPAxis::stLinear);
        addGraph(composite_plots[i], p[0].times, plot_vals[i]);
    }

    std::cout<<"Obliczono"<<std::endl;
}

void MainWindow::on_clearButton_clicked()
{

    for(auto&& p: plots){
        p->clearGraphs();
        p->replot();
    }
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
    ui->hexRadioButton->setChecked(false);
    ui->cubicRadioButton->setChecked(true);
}

void MainWindow::on_hexRadioButton_clicked(){
    ui->hexRadioButton->setChecked(true);
    ui->cubicRadioButton->setChecked(false);
}
