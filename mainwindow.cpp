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

    std::vector<QCustomPlot*> plots{ui->kappaPlot, ui->lambdaPlot, ui->etaPlot,
                                    ui->lambda2Plot, ui->lambda3Plot, ui->lambda4Plot};

    QVector<double> x=QVector<double>::fromStdVector(p.times);

    for(size_t i=0; i<plots.size(); i++){
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
    }
}

PlotSet MainWindow::integrate(){
    std::vector<double> coords({ui->kappaBox->value(), ui->lambdaBox->value(), ui->etaBox->value(),
                                ui->lambda2Box->value(), ui->lambda3Box->value(), ui->lambda4Box->value()});
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
    PlotSet p = this->integrate();
    std::ofstream file(ui->filenameEdit->text().toStdString());
    file << p;
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

void MainWindow::on_lambda2Box_valueChanged(double)
{
    this->update();
}

void MainWindow::on_lambda3Box_valueChanged(double)
{
    this->update();
}

void MainWindow::on_lambda4Box_valueChanged(double)
{
    this->update();
}

