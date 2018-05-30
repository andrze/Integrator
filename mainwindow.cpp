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
#include "physics.cpp"
#include <iostream>

extern const EquationSet equations;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    integrator = Integrator(equations);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::update(){
    PlotSet p = this->integrate();

    std::vector<std::vector<double> > results = p.transpose();

    std::vector<QCustomPlot*> plots{ui->kappaPlot, ui->uPlot, ui->lambdaPlot, ui->etaPlot,
                                    ui->kappa2Plot, ui->u2Plot, ui->lambda2Plot, ui->eta2Plot};

    QVector<double> x=QVector<double>::fromStdVector(p.times);

    std::vector<double> res;

    for(size_t i=0; i<plots.size()/2; i++){
        auto plot = plots[i];
        std::vector<double> vals = results[i];
        QVector<double> y = QVector<double>::fromStdVector(vals);
        plot->addGraph();
        plot->graph(0)->setData(x, y);
        plot->xAxis->setRange(x.first(),x.last());
        auto min = std::min_element(vals.begin(), vals.end());
        auto max = std::max_element(vals.begin(), vals.end());
        plot->yAxis->setRange((*min)-0.25, (*max)+0.25),
        plot->replot();

        //Unscaled plot
        plot = plots[i+plots.size()/2];
        for(int j=0; j<x.size(); j++){
            vals[j] *= exp(integrator.equations.scale[i] * x[j]);
        }
        y = QVector<double>::fromStdVector(vals);
        plot->addGraph();
        plot->graph(0)->setData(x, y);
        plot->xAxis->setRange(x.first(),x.last());
        min = std::min_element(vals.begin(), vals.end());
        max = std::max_element(vals.begin(), vals.end());
        plot->yAxis->setRange((*min)-0.25, (*max)+0.25),
        plot->replot();

        res.push_back(results[i][0]);
        res.push_back(vals[vals.size()-1]);
    }

    this->results.vals.push_back(RealVector(res));
}

PlotSet MainWindow::integrate(){
    std::vector<double> coords({ui->kappaBox->value(), ui->uBox->value(),
                                ui->lambdaBox->value(), ui->etaBox->value()});
    RealVector start(coords);

    return integrator.integrate(ui->startTimeBox->value(),
                                ui->endTimeBox->value(),
                                ui->deltaTBox->value(), start);
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

void MainWindow::on_kappaBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_lambdaBox_valueChanged(double)
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

void MainWindow::on_etaBox_valueChanged(double)
{
    this->update();
}


void MainWindow::on_uBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_ksenseBox_valueChanged(int arg1)
{
    this->ui->kappaBox->setSingleStep(pow(10.,-arg1));
}

void MainWindow::on_usenseBox_valueChanged(int arg1)
{
    this->ui->uBox->setSingleStep(pow(10.,-arg1));
}



void MainWindow::on_lsenseBox_valueChanged(int arg1)
{
    this->ui->lambdaBox->setSingleStep(pow(10.,-arg1));
}
