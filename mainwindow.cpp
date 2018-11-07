#include <vector>
#include <functional>
#include <fstream>
#include <QString>
#include <QFileDialog>
#include <algorithm>
#include <cmath>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "integrator.h"
#include "realvector.h"
#include "physics.h"
#include <iostream>

void setupPlot(QCustomPlot* plot, QVector<double> x, std::vector<double> vals){
    QVector<double> y = QVector<double>::fromStdVector(vals);
    plot->addGraph();
    plot->graph(0)->setData(x, y);
    plot->xAxis->setRange(x.first(),x.last());
    auto min = *std::min_element(vals.begin(), vals.end());
    auto max = *std::max_element(vals.begin(), vals.end());
    double diff = max - min;
    if(diff == 0){
        diff = 1.;
    }
    plot->yAxis->setRange(min-diff/4, max+diff/4);
    plot->replot();
}


void setupPlot(QCustomPlot* plot, QVector<double> x, std::vector<double> vals, std::vector<double> vals2){
    QVector<double> y = QVector<double>::fromStdVector(vals);
    plot->addGraph();
    plot->graph(0)->setPen(QPen(Qt::blue));
    plot->graph(0)->setData(x, y);
    plot->xAxis->setRange(x.first(),x.last());
    auto min1 = std::min_element(vals.begin(), vals.end());
    auto max1 = std::max_element(vals.begin(), vals.end());

    QVector<double> y2 = QVector<double>::fromStdVector(vals2);
    plot->addGraph();
    plot->graph(1)->setPen(QPen(Qt::red));
    plot->graph(1)->setData(x, y2);

    auto min2 = std::min_element(vals2.begin(), vals2.end());
    auto max2 = std::max_element(vals2.begin(), vals2.end());
    double min = std::min(*min1, *min2);
    double max = std::max(*max1, *max2);
    double diff = max - min;
    if(diff == 0){
        diff = 1.;
    }
    plot->yAxis->setRange(min-diff/4, max+diff/4),
    plot->replot();
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    EquationSet equations;
    integrator = Integrator(equations);
    ui->asenseBox->valueChanged(ui->asenseBox->value());
    ui->usenseBox->valueChanged(ui->usenseBox->value());
    ui->tsenseBox->valueChanged(ui->tsenseBox->value());
    this->update();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::update(){
    double a = std::pow(ui->alphaBox->value(),2);
    double u = ui->uBox->value();
    double tau = ui->tauBox->value();
    double Z = ui->ZSpinBox->value();
    double Y = ui->YSpinBox->value();

    ui->msBox->setValue(a*u);
    ui->mpBox->setValue(a*tau);
    ui->zsSpinBox->setValue(a*Y+Z);
    ui->zpSpinBox->setValue(Z);

}

PlotSet MainWindow::integrate(){
    std::vector<double> coords({ui->alphaBox->value(),
                                ui->msBox->value(), ui->mpBox->value(),
                                ui->zsSpinBox->value(), ui->zpSpinBox->value(),
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
    file << results;
}

void MainWindow::on_alphaBox_valueChanged(double)
{
    this->update();
}


void MainWindow::on_uBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_tauBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_ZSpinBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_YSpinBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_startTimeBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_endTimeBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_deltaTBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_asenseBox_valueChanged(int arg1)
{
    this->ui->alphaBox->setSingleStep(pow(10.,-arg1));
}

void MainWindow::on_usenseBox_valueChanged(int arg1)
{
    this->ui->uBox->setSingleStep(pow(10.,-arg1));
}



void MainWindow::on_tsenseBox_valueChanged(int arg1)
{
    this->ui->tauBox->setSingleStep(pow(10.,-arg1));
}

/*
void MainWindow::on_fpButton_clicked()
{
    double lower = 0., current = ui->kappaBox->value(), upper;
    int steps_left = 64;

    std::vector<double> coords({current, ui->uBox->value(),
                                ui->lambdaBox->value(), 1});
    RealVector point(coords);
    Integrator integ(this->integrator);
    PlotSet res = integ.integrate(0, -20, -0.01,
                                  point, ui->dSpinBox->value());
    RealVector end = res.vals.back();
    RealVector der = integ.equations.evaluate(end,
                                              ui->dSpinBox->value());

    while(der[0]>=0.){
        current = 2*(std::abs(current) + 0.01);
        point[0] = current;

        res = integ.integrate(0, -20, -0.01,
                              point, ui->dSpinBox->value());

        steps_left--;

        end = res.vals.back();
        der = integ.equations.evaluate(end,
                                       ui->dSpinBox->value());
        if(steps_left==0){
            ui->kappaBox->setValue(current);
            return;
        }
    }

    upper = current;

    for(;steps_left>0; steps_left--){
        if(der[0] < 0.){

            upper = current;
            current = (lower+current)/2;
        } else {
            lower = current;
            current = (upper+current)/2;
        }
        point[0] = current;
        res = integ.integrate(0, -20, -0.01,
                              point, ui->dSpinBox->value());

        end = res.vals.back();
        der = integ.equations.evaluate(end,
                                       ui->dSpinBox->value());
    }

    ui->kappaBox->setValue(current);
    return;
}*/

void MainWindow::on_pushButton_clicked(){


    PlotSet p = this->integrate();
    std::vector<std::vector<double> > results = p.transpose();
    std::array<std::vector<double>, 6> plot_vals;

    std::vector<QCustomPlot*> autoplots{ui->alphaPlot, ui->msPlot, ui->mpPlot};
    std::vector<QCustomPlot*> composite_plots{ui->uPlot, ui->tauPlot, ui->yPlot, ui->stiffPlot, ui->etaPlot, ui->universalPlot};

    QVector<double> x;
    for(auto&& t: p.times){
        x.push_back(-log(t));
    }


    for(size_t i=0; i<autoplots.size(); i++){
        autoplots[i]->yAxis->setScaleType(QCPAxis::stLogarithmic);
        setupPlot(autoplots[i], x, results[i]);
    }
    ui->zPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    setupPlot(ui->zPlot, x, results[3], results[4]);

    for(int i=0; i<x.size(); i++){
        double a2 = std::pow(results[0][i],2);
        plot_vals[0].push_back((results[1][i])/a2);
        plot_vals[1].push_back((results[2][i])/a2);
        plot_vals[2].push_back((results[3][i]-results[4][i])/a2);
        plot_vals[3].push_back(a2*results[4][i]);
        plot_vals[4].push_back(-p.times[i]*p.derivatives[i][4]/results[4][i]);
        plot_vals[5].push_back(-p.times[i]*p.derivatives[i][4]*a2/results[5][i]);
    }

    for(size_t i=0; i<composite_plots.size(); i++){
        composite_plots[i]->yAxis->setScaleType(QCPAxis::stLinear);
        setupPlot(composite_plots[i], x, plot_vals[i]);
    }

    std::cout<<"Obliczono"<<std::endl;
}
